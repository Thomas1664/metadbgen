import sqlite3
import asyncio
import json
import threading
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Optional, Type
from pydantic import BaseModel


class AsyncSQLiteCache:
    _lock = threading.Lock()

    def __init__(self, db_path: str = "cache.sqlite"):
        self.db_path = Path(db_path)
        self._init_db()

    def _init_db(self) -> None:
        with self._lock, sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS cache (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                )
                """
            )

    def _to_primitive(self, obj: Any) -> Any:
        """
        Convert obj into JSON-serializable primitives.
        Crucial: for BaseModel instances, use model_dump(by_alias=True).
        This recursively converts lists/dicts/sets/tuples too.
        """
        # Pydantic model -> use alias-aware dump
        if isinstance(obj, BaseModel):
            # model_dump(by_alias=True) returns dicts/lists using validation/serialization aliases
            return obj.model_dump(by_alias=True)

        # dict: convert values
        if isinstance(obj, dict):
            return {k: self._to_primitive(v) for k, v in obj.items()}

        # iterable containers -> convert items (preserve as list)
        if isinstance(obj, (list, tuple, set)):
            return [self._to_primitive(v) for v in obj]

        # primitive / json-serializable object -> return as-is
        return obj

    def _json_dumps(self, obj: Any) -> str:
        prim = self._to_primitive(obj)
        # Use separators to make stored JSON compact and ensure_ascii=False to keep unicode
        return json.dumps(prim, ensure_ascii=False, separators=(",", ":"))

    def get(self, key: str) -> Optional[Any]:
        with self._lock, sqlite3.connect(self.db_path) as conn:
            row = conn.execute("SELECT value FROM cache WHERE key=?", (key,)).fetchone()
            return json.loads(row[0]) if row else None

    def set(self, key: str, value: Any) -> None:
        json_val = self._json_dumps(value)
        with self._lock, sqlite3.connect(self.db_path) as conn:
            conn.execute("REPLACE INTO cache (key, value) VALUES (?, ?)", (key, json_val))
            conn.commit()

    def decorator(self, key_func: Callable[..., str], model_cls: Optional[Type[BaseModel]] = None):
        """
        key_func(*args, **kwargs) -> str: stable cache key (string)
        model_cls: the element model class if the cached value is a list of models
        """
        def wrapper(fn):
            @wraps(fn)
            async def async_wrapped(*args, **kwargs):
                key = key_func(*args, **kwargs)
                cached = self.get(key)
                if cached is not None:
                    if model_cls:
                        # cached is primitive alias-based JSON -> reconstruct back to model(s)
                        if isinstance(cached, list):
                            return [model_cls.model_validate(item) for item in cached]
                        if isinstance(cached, dict):
                            return model_cls.model_validate(cached)
                    return cached

                # Cache miss -> call the async function
                result = await fn(*args, **kwargs)
                if result is None:
                    return None

                # Persist result in background thread so we don't block the event loop
                loop = asyncio.get_running_loop()
                await loop.run_in_executor(None, lambda: self.set(key, result))
                return result

            return async_wrapped
        return wrapper
