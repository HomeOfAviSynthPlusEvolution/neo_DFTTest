from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
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


def render_vs_script(plugin_path: str, source: dict, params: dict) -> str:
    if source["type"] != "blank":
        raise ValueError(f"unsupported source type: {source['type']}")

    param_text = ", ".join(
        f"{key}={_format_vs_value(params[key])}" for key in sorted(params)
    )
    if param_text:
        param_text = ", " + param_text

    return "\n".join(
        [
            "import vapoursynth as vs",
            "core = vs.core",
            f"core.std.LoadPlugin(path={json.dumps(plugin_path)})",
            (
                "src = core.std.BlankClip("
                f"width={source['width']}, "
                f"height={source['height']}, "
                f"length={source['length']}, "
                f"format=vs.{source['format']}, "
                f"color={_format_vs_value(source['color'])})"
            ),
            f"clip = core.neo_dfttest.DFTTest(src{param_text})",
            "clip.set_output()",
            "",
        ]
    )


def render_avs_script(plugin_path: str, source: dict, params: dict) -> str:
    if source["type"] != "blank":
        raise ValueError(f"unsupported source type: {source['type']}")

    param_text = ", ".join(
        f"{key}={_format_avs_value(params[key])}" for key in sorted(params)
    )
    if param_text:
        param_text = ", " + param_text

    return "\n".join(
        [
            f"LoadPlugin({json.dumps(plugin_path)})",
            (
                "src = BlankClip("
                f"width={source['width']}, "
                f"height={source['height']}, "
                f"length={source['length']}, "
                f"pixel_type={json.dumps(source['format'])}, "
                f"color={_format_avs_value(source['color'])})"
            ),
            f"return neo_dfttest(src{param_text})",
            "",
        ]
    )


def _format_vs_value(value) -> str:
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, list):
        return "[" + ", ".join(_format_vs_value(item) for item in value) + "]"
    return repr(value)


def _format_avs_value(value) -> str:
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        return json.dumps(", ".join(str(item) for item in value))
    return repr(value)
