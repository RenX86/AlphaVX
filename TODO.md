1. ~~**Concurrency / Async Batching**: Use `asyncio` or ThreadPoolExecutor in `score_batch` to score multiple variants concurrently and speed up batch jobs.~~
2. ~~**Annotated VCF Export**: Add support for outputting an annotated `.vcf` file with scores embedded in the `INFO` column to integrate with bioinformatics pipelines.~~
3. **Local Pre-filtering**: Add options to filter variants locally by BED file (regions) or allele frequency (gnomAD) before querying the API, to save time/quota.
4. ~~**Interactive HTML Report**: Upgrade the static HTML report to use DataTables.js and interactive Plotly graphs.~~
5. ~~**CI/CD Pipeline**: Add GitHub Actions workflows for automated testing, linting, and PyPI publishing.~~