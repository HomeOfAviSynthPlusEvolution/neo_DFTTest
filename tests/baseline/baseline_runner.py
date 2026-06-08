from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Iterable


@dataclass(frozen=True)
class PlaneBytes:
    name: str
    width: int
    height: int
    bytes_per_sample: int
    stride: int
    data: bytes

    @property
    def row_bytes(self) -> int:
        return self.width * self.bytes_per_sample


def hash_planes(planes: Iterable[PlaneBytes]) -> str:
    digest = hashlib.sha256()
    for plane in planes:
        row_bytes = plane.row_bytes
        for y in range(plane.height):
            start = y * plane.stride
            digest.update(plane.data[start:start + row_bytes])
    return digest.hexdigest()


def select_cases(cases: list[dict], tier: str) -> list[dict]:
    if tier == "all":
        return list(cases)
    return [case for case in cases if case.get("tier") == tier]
