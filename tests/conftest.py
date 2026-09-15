"""Shared test fixtures for the AlphaVX test suite.

Creates mock modules for the alphagenome SDK when it's not installed,
so scorer tests can run without the real dependency.
"""

from __future__ import annotations

import sys
from unittest.mock import MagicMock

import pytest


def _ensure_alphagenome_mocks():
    """Install mock alphagenome modules in sys.modules if not already present."""
    modules_to_mock = [
        "alphagenome",
        "alphagenome.data",
        "alphagenome.data.genome",
        "alphagenome.models",
        "alphagenome.models.dna_client",
        "alphagenome.models.variant_scorers",
    ]
    for mod_name in modules_to_mock:
        if mod_name not in sys.modules:
            mock_mod = MagicMock()
            if mod_name == "alphagenome.models.variant_scorers":
                mock_mod.RECOMMENDED_VARIANT_SCORERS = {}
            sys.modules[mod_name] = mock_mod
    
    # Link submodules to parents so 'from alphagenome.models import variant_scorers' works
    if "alphagenome" in sys.modules and isinstance(sys.modules["alphagenome"], MagicMock):
        sys.modules["alphagenome"].data = sys.modules.get("alphagenome.data")
        sys.modules["alphagenome"].models = sys.modules.get("alphagenome.models")
        
    if "alphagenome.data" in sys.modules and isinstance(sys.modules["alphagenome.data"], MagicMock):
        sys.modules["alphagenome.data"].genome = sys.modules.get("alphagenome.data.genome")
        
    if "alphagenome.models" in sys.modules and isinstance(sys.modules["alphagenome.models"], MagicMock):
        sys.modules["alphagenome.models"].dna_client = sys.modules.get("alphagenome.models.dna_client")
        sys.modules["alphagenome.models"].variant_scorers = sys.modules.get("alphagenome.models.variant_scorers")

# Install mocks at import time so scorer.py can be imported
_ensure_alphagenome_mocks()
