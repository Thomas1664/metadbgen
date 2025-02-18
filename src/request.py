from concurrent.futures import ThreadPoolExecutor
import httpx
import os
from dotenv import load_dotenv
import asyncio
from bs4 import BeautifulSoup
from pathlib import Path

class HttpClient:
    def __init__(self) -> None:
        load_dotenv()
        self._base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self._client = httpx.AsyncClient(params={
            'api_key': os.getenv('NCBI_API_KEY'),
            'email': os.getenv('NCBI_API_EMAIL'),
            'tool': os.getenv('NCBI_API_TOOL')
        })

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
        return (query_key, web_env, count)
    
    async def efetch(self, query_key: str, web_env: str, ret_start: str, out_file: str, db = 'nuccore'):
        response = await self._client.get(f'{self._base_url}efetch.fcgi', params={
            'db': db,
            'query_key': query_key,
            'WebEnv': web_env,
            'retmode': 'text',
            'rettype': 'gb',
            'retmax': '500',
            'retstart': ret_start
        })
        response.raise_for_status()

        path = Path(out_file)
        with open(f'{path.stem}_{ret_start}{path.suffix}', 'w') as file:
            print(response.text, file=file)
    
    async def fetch_results(self, query_key: str, web_env: str, n_items: str, out_file, db = 'nuccore'):
        retstart_values = list(range(0, int(n_items), 500))
        #with ThreadPoolExecutor(max_workers=5) as executor:
        for start in retstart_values:
            #executor.map(lambda start: asyncio.run(self.efetch(query_key, web_env, start, out_file, db)), retstart_values)
            await self.efetch(query_key, web_env, str(start), out_file, db)

async def main():
    client = HttpClient()
    response = await client.esearch('cox2 AND oomycota[organism]')
    
    await client.fetch_results(*response, 'data/cox2.gb')
    
asyncio.run(main())
