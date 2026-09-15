"""Tests for the configuration module.

Covers:
- YAML loading with all fields
- Environment variable resolution
- key_env YAML field support (bug fix verification)
- Default values
- Missing API key error
"""

from __future__ import annotations

import os
from pathlib import Path
from textwrap import dedent

import pytest

from alphavx.config import Config, load_config, DEFAULT_MODALITIES


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def yaml_file(tmp_path: Path) -> Path:
    """Create a minimal valid alphavx.yaml."""
    content = dedent("""\
        api:
          key_env: ALPHAVX_API_KEY
          max_retries: 5
          retry_delay: 10.0

        scoring:
          modalities:
            - RNA_SEQ
            - DNASE
          quantile_threshold: 0.999
          sequence_length: 524288

        output:
          cache: false
    """)
    path = tmp_path / "alphavx.yaml"
    path.write_text(content)
    return path


# ---------------------------------------------------------------------------
# Config dataclass defaults
# ---------------------------------------------------------------------------

class TestConfigDefaults:
    """Verify Config dataclass defaults."""

    def test_default_modalities(self):
        c = Config()
        assert c.modalities == list(DEFAULT_MODALITIES)

    def test_default_values(self):
        c = Config()
        assert c.max_retries == 3
        assert c.retry_delay == 5.0
        assert c.sequence_length == 2**20
        assert c.quantile_threshold == 0.995
        assert c.cache_enabled is True
        assert c.api_key == ""


# ---------------------------------------------------------------------------
# YAML loading
# ---------------------------------------------------------------------------

class TestYamlLoading:
    """Tests for load_config with YAML files."""

    def test_loads_all_yaml_fields(self, yaml_file: Path, monkeypatch):
        monkeypatch.setenv("ALPHAVX_API_KEY", "test-key-123")

        config = load_config(yaml_file)

        assert config.max_retries == 5
        assert config.retry_delay == 10.0
        assert config.modalities == ["RNA_SEQ", "DNASE"]
        assert config.quantile_threshold == 0.999
        assert config.sequence_length == 524288
        assert config.cache_enabled is False
        assert config.api_key == "test-key-123"

    def test_partial_yaml(self, tmp_path: Path, monkeypatch):
        """YAML with only some fields should leave others at defaults."""
        monkeypatch.setenv("ALPHAVX_API_KEY", "key")
        yaml = tmp_path / "partial.yaml"
        yaml.write_text("api:\n  max_retries: 10\n")

        config = load_config(yaml)

        assert config.max_retries == 10
        assert config.retry_delay == 5.0  # default
        assert config.modalities == list(DEFAULT_MODALITIES)  # default

    def test_empty_yaml(self, tmp_path: Path, monkeypatch):
        """Empty YAML should produce default config."""
        monkeypatch.setenv("ALPHAVX_API_KEY", "key")
        yaml = tmp_path / "empty.yaml"
        yaml.write_text("")

        config = load_config(yaml)

        assert config.max_retries == 3
        assert config.modalities == list(DEFAULT_MODALITIES)


# ---------------------------------------------------------------------------
# API key resolution
# ---------------------------------------------------------------------------

class TestApiKeyResolution:
    """Tests for API key lookup from environment variables."""

    def test_alphavx_key_preferred(self, monkeypatch):
        monkeypatch.setenv("ALPHAVX_API_KEY", "alphavx-key")
        monkeypatch.setenv("ALPHAGENOME_API_KEY", "alphagenome-key")

        config = load_config(None)

        assert config.api_key == "alphavx-key"

    def test_alphagenome_fallback(self, monkeypatch):
        monkeypatch.delenv("ALPHAVX_API_KEY", raising=False)
        monkeypatch.setenv("ALPHAGENOME_API_KEY", "fallback-key")

        config = load_config(None)

        assert config.api_key == "fallback-key"

    def test_no_key_raises(self, monkeypatch):
        monkeypatch.delenv("ALPHAVX_API_KEY", raising=False)
        monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)

        with pytest.raises(ValueError, match="No API key found"):
            load_config(None)

    def test_key_env_yaml_field_respected(self, tmp_path: Path, monkeypatch):
        """key_env in YAML should be used to look up the API key.

        This verifies the fix for the decorative key_env bug.
        """
        monkeypatch.delenv("ALPHAVX_API_KEY", raising=False)
        monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
        monkeypatch.setenv("MY_CUSTOM_KEY_VAR", "custom-key-value")

        yaml = tmp_path / "custom.yaml"
        yaml.write_text("api:\n  key_env: MY_CUSTOM_KEY_VAR\n")

        config = load_config(yaml)

        assert config.api_key == "custom-key-value"

    def test_key_env_falls_back_to_defaults(self, tmp_path: Path, monkeypatch):
        """If key_env var is empty/unset, should still try ALPHAVX/ALPHAGENOME."""
        monkeypatch.delenv("MISSING_VAR", raising=False)
        monkeypatch.setenv("ALPHAVX_API_KEY", "default-key")

        yaml = tmp_path / "custom.yaml"
        yaml.write_text("api:\n  key_env: MISSING_VAR\n")

        config = load_config(yaml)

        assert config.api_key == "default-key"


# ---------------------------------------------------------------------------
# Config.scoring_fingerprint (cache invalidation fix)
# ---------------------------------------------------------------------------

class TestScoringFingerprint:
    """Tests for the config fingerprint used in cache invalidation."""

    def test_fingerprint_exists(self):
        """Config should have a scoring_fingerprint property."""
        c = Config(api_key="k", modalities=["RNA_SEQ"])
        assert hasattr(c, "scoring_fingerprint")

    def test_same_config_same_fingerprint(self):
        """Identical scoring config should produce identical fingerprints."""
        c1 = Config(api_key="k1", modalities=["RNA_SEQ", "DNASE"], sequence_length=1024)
        c2 = Config(api_key="k2", modalities=["RNA_SEQ", "DNASE"], sequence_length=1024)
        assert c1.scoring_fingerprint == c2.scoring_fingerprint

    def test_different_modalities_different_fingerprint(self):
        c1 = Config(api_key="k", modalities=["RNA_SEQ"])
        c2 = Config(api_key="k", modalities=["RNA_SEQ", "DNASE"])
        assert c1.scoring_fingerprint != c2.scoring_fingerprint

    def test_different_sequence_length_different_fingerprint(self):
        c1 = Config(api_key="k", sequence_length=1024)
        c2 = Config(api_key="k", sequence_length=2048)
        assert c1.scoring_fingerprint != c2.scoring_fingerprint

    def test_modality_order_irrelevant(self):
        """Modality order shouldn't affect fingerprint (sorted internally)."""
        c1 = Config(api_key="k", modalities=["DNASE", "RNA_SEQ"])
        c2 = Config(api_key="k", modalities=["RNA_SEQ", "DNASE"])
        assert c1.scoring_fingerprint == c2.scoring_fingerprint

    def test_api_key_doesnt_affect_fingerprint(self):
        """API key is not a scoring parameter — shouldn't change fingerprint."""
        c1 = Config(api_key="key-1", modalities=["RNA_SEQ"])
        c2 = Config(api_key="key-2", modalities=["RNA_SEQ"])
        assert c1.scoring_fingerprint == c2.scoring_fingerprint
