from dataclasses import dataclass
import httpx
import os
from dotenv import load_dotenv
from bs4 import BeautifulSoup
from pathlib import Path


@dataclass
class ESearchResult:
    query_key: str
    web_env: str
    result_count: int

class HttpClient:
    def __init__(self) -> None:
        load_dotenv()
        self._base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self._client = httpx.AsyncClient(params={
            'api_key': os.getenv('NCBI_API_KEY'),
            'email': os.getenv('NCBI_API_EMAIL'),
            'tool': os.getenv('NCBI_API_TOOL')
        }, timeout=httpx.Timeout(30, read=30))

    async def esearch(self, query: str, db='nuccore'):
        response = await self._client.get(f'{self._base_url}esearch.fcgi', params={
            'db': db,
            'usehistory': 'y',
            'term': query
        })
        response.raise_for_status()
        bs = BeautifulSoup(response.text, features="xml")
        print(response.text)
        query_key: str = bs.find('QueryKey').string
        web_env: str = bs.find('WebEnv').string
        count: str = bs.find('Count').string
        return ESearchResult(query_key, web_env, int(count))

    async def efetch_all(self, esearch_result: ESearchResult, out_file: str, db='nuccore', rettype: str='gb'):
        params = params={
                'db': db,
                'query_key': esearch_result.query_key,
                'WebEnv': esearch_result.web_env,
                'retmode': 'text',
                'rettype': rettype
            }
        async with self._client.stream("GET", f'{self._base_url}efetch.fcgi', params=params) as response:
            response.raise_for_status()
            path = Path(out_file)
            out_file = f'{path.parent}/{path.stem}{path.suffix}'
            with open(out_file, "wb") as f:
                async for chunk in response.aiter_bytes():
                    f.write(chunk)    

async def fetch_sequences(query_str: str, file_prefix: str):
    client = HttpClient()
    response = await client.esearch(query_str)
    await client.efetch_all(response, file_prefix)
