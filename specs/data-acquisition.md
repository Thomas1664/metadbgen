# Data Acquisition

## Sequences

Sequences are downloaded using the [NCBI ENTREZ API](https://www.ncbi.nlm.nih.gov/books/NBK25501/). It allows automated [download of large datasets](https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large).

## Sequence metadata

### Taxonomical data

Previously, taxonomical data were obtained from [Index Fungorum](https://www.indexfungorum.org). There is no API, so we would to scrape the data from their website. Furthermore, There doesn't seem to be information about the full lineage available right away.

Therefore, we use the [taxonomy browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) from NCBI. It offers an API and lineage information.

The taxonomical information will be used to validate the organism assignments of downloaded sequences.
