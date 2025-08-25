# Downloading sequence data

File: [src/request.py]()

## HttpClient

The class **HttpClient** handles all web requests to the [NCBI entrez API](https://www.ncbi.nlm.nih.gov/books/NBK25501/) using the asynchronous HTTP client from [`httpx`](https://www.python-httpx.org/).

All requests to the entrez API require 3 URL parameters: `api_key`, `email`, and `tool`. The values of these parameters are determined by the environment variables `NCBI_API_KEY`, `NCBI_API_EMAIL`, and `NCBI_API_TOOL`, respectively. An API key is required and increases the rate limit of the entrez API from 3 to 10 requests per second. You can obtain an API key using your NCBI account. It is the responsibility of the user to ensure the rate limit is not exceeded. `NCBI_API_EMAIL` should be set to your personal email address. `NCBI_API_TOOL` should be set to `cox_auto_downloader`.

### esearch

`esearch` is used to find sequences in the NCBI databases based on a query. The query is specified by setting the URL parameter `term` to the query string. Spechial characters like spaces or brackets don't need to be escaped. The nucleotide database is used by default setting the `db` URL parameter to `nuccore`. `usehistory` is set to `y` to allow use of the history server to fetch the results.

A successful `esearch` request returns a XML response containing a `QueryKey` and a `WebEnv` which can be used to fetch the results from the history server using `efetch`. Finally, the hit count is returned in `Count`

### efetch

`efetch` is used to retrieve the results from the history server using `query_key` and `WebEnv` which were returned by a previous `esearch` request. The results are obtained as GenBank records by setting the URL parameters `retmode` to `text` and `rettype` to `gb`.

## Fetching Cox2 sequences

All available Cox2 sequences in Oomycetes were obtained from the NCBI nucleotide database using the query `cox2 AND oomycota[organism]`. However, data analysis revealed different spellings of Cox2. Therefore, the alternative spellings of Cox2 should also be considered in the search query.
