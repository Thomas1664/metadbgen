from collections.abc import Callable
from dataclasses import dataclass
import httpx
from httpx_retries import RetryTransport, Retry
import os
from dotenv import load_dotenv
from bs4 import BeautifulSoup
from pathlib import Path
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential
import asyncio


@dataclass
class ESearchResult:
    query_key: str
    web_env: str
    result_count: int
    text_response: str

class HttpClient:
    def __init__(self) -> None:
        load_dotenv()
        api_key = os.getenv('NCBI_API_KEY')
        if api_key is None:
            raise RuntimeError("Missing NCBI API key!")
        self._base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self._client = httpx.AsyncClient(
            params={
                'api_key': api_key,
                'email': os.getenv('NCBI_API_EMAIL'),
                'tool': os.getenv('NCBI_API_TOOL')
            },
            timeout=httpx.Timeout(60, read=60),
            transport=RetryTransport(retry=Retry(total=5, backoff_factor=0.5, status_forcelist=[400, 429, 502, 503, 504]))
        )

    @retry(
        retry=retry_if_exception_type(RuntimeError),
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=10),
        reraise=True
    )
    async def esearch(self, query: str, db='nuccore'):
        response = await self._client.get(f'{self._base_url}esearch.fcgi', params={
            'db': db,
            'usehistory': 'y',
            'term': query
        })
        response.raise_for_status()
        bs = BeautifulSoup(response.text, features="xml")
        maybe_query_key = bs.find('QueryKey')
        maybe_web_env = bs.find('WebEnv')
        maybe_count = bs.find('Count')
        if maybe_query_key is None or maybe_web_env is None or maybe_count is None:
            raise RuntimeError(f"Failed to parse esearch result:\n{response.text}")
        return ESearchResult(maybe_query_key.string, maybe_web_env.string, int(maybe_count.string), response.text)

    async def efetch_all(self, esearch_result: ESearchResult, out_file: str | Path, db='nuccore', rettype: str='gb'):
        out_files: set[str] = set()
        # Can only fetch 10K entries per request!
        MAX_ENTRIES = 10000
        path = Path(out_file)
        for start in range(0, esearch_result.result_count, MAX_ENTRIES):
            print(f'Downloading {start} of {esearch_result.result_count} records...')
            params = {
                'db': db,
                'query_key': esearch_result.query_key,
                'WebEnv': esearch_result.web_env,
                'retmode': 'text',
                'rettype': rettype,
                'retstart': start,
                'retmax': MAX_ENTRIES
            }
            async with self._client.stream("GET", f'{self._base_url}efetch.fcgi', params=params) as response:
                response.raise_for_status()
                current_out_file = f'{path.parent}/{path.stem}_{start}{path.suffix}'
                with open(current_out_file, "wb") as f:
                    async for chunk in response.aiter_bytes():
                        f.write(chunk)
                out_files.add(current_out_file)
        # Merge GenBank files back into one file.
        # This is required because output file is tracked by snakemake.
        if rettype == 'gb':
            with open(out_file, 'w') as single_output:
                for file in out_files:
                    with open(file, 'r') as part_f:
                        single_output.write(part_f.read())
            return out_file
        return out_files

    # ids: list of ids obtained from esearch; make sure there are no None values
    # out_file_pattern: format string that allows insertion of a integer representing each batch
    async def efetch_by_ids(self, db: str, ids: list[str], callback: Callable, batch_size=100, ret_type='xml') -> list[str]:
        outputs = []
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i+batch_size]
            response = await self._client.get(f'{self._base_url}efetch.fcgi', params={
                'db': db,
                "id": ",".join(batch),
                "retmode": "text",
                "rettype": ret_type,
            })
            response.raise_for_status()
            outputs.append(callback(response.text))
            # rate limit: 10 requests per second
            await asyncio.sleep(0.1)
        return outputs
    
    async def efetch_by_id(self, db: str, id: str, ret_mode='xml') -> str:
        response = await self._client.get(
            f'{self._base_url}efetch.fcgi',
            params={
                'db': db,
                "id": id,
                "retmode": ret_mode,
                "rettype": "text"
            }
        )
        response.raise_for_status()
        return response.text

async def fetch_sequences(query_str: str, file_prefix: str):
    client = HttpClient()
    response = await client.esearch(query_str)
    await client.efetch_all(response, file_prefix)

class NcbiApiClient:
    def __init__(self) -> None:
        load_dotenv()
        api_key = os.getenv('NCBI_API_KEY')
        if api_key is None:
            raise RuntimeError("Missing NCBI API key!")
        self._base_url = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
        self._client = httpx.AsyncClient(
            headers={
                'api-key': api_key,
            },
            timeout=httpx.Timeout(60, read=60),
            transport=RetryTransport(retry=Retry(total=5, backoff_factor=0.5))
        )
    
    async def post(self, dataset: str, endpoint: str, json: dict):
        response = await self._client.post(f'{self._base_url}/{dataset}/{endpoint}', json=json)
        response.raise_for_status()
        return response.json()
