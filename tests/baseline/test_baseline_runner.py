import hashlib
import unittest

from tests.baseline import baseline_runner


class CanonicalHashTests(unittest.TestCase):
    def test_hashes_only_valid_row_bytes(self):
        planes = [
            baseline_runner.PlaneBytes(
                name="Y",
                width=3,
                height=2,
                bytes_per_sample=1,
                stride=5,
                data=b"abcXXdefYY",
            )
        ]

        result = baseline_runner.hash_planes(planes)

        self.assertEqual(result, hashlib.sha256(b"abcdef").hexdigest())

    def test_hashes_planes_in_supplied_order(self):
        planes = [
            baseline_runner.PlaneBytes("Y", 2, 1, 1, 2, b"ab"),
            baseline_runner.PlaneBytes("U", 2, 1, 1, 2, b"cd"),
            baseline_runner.PlaneBytes("V", 2, 1, 1, 2, b"ef"),
        ]

        result = baseline_runner.hash_planes(planes)

        self.assertEqual(result, hashlib.sha256(b"abcdef").hexdigest())


class CaseSelectionTests(unittest.TestCase):
    def test_filters_cases_by_tier(self):
        cases = [
            {"id": "smoke_default", "tier": "smoke"},
            {"id": "compat_yuv", "tier": "compat"},
            {"id": "full_hd", "tier": "full"},
        ]

        selected = baseline_runner.select_cases(cases, "smoke")

        self.assertEqual([case["id"] for case in selected], ["smoke_default"])

    def test_all_tier_keeps_all_cases(self):
        cases = [
            {"id": "smoke_default", "tier": "smoke"},
            {"id": "compat_yuv", "tier": "compat"},
        ]

        selected = baseline_runner.select_cases(cases, "all")

        self.assertEqual([case["id"] for case in selected], ["smoke_default", "compat_yuv"])


class ScriptRenderingTests(unittest.TestCase):
    def test_renders_vs_filter_call_with_arrays_and_strings(self):
        source = {
            "type": "blank",
            "format": "GRAY8",
            "width": 64,
            "height": 48,
            "length": 9,
            "color": [64],
        }
        params = {"tbsize": 1, "sigma": 8.0, "planes": [0], "mode": "debug"}

        script = baseline_runner.render_vs_script("/tmp/plugin.so", source, params)

        self.assertIn('core.std.LoadPlugin(path="/tmp/plugin.so")', script)
        self.assertIn("core.std.BlankClip(width=64, height=48, length=9, format=vs.GRAY8, color=[64])", script)
        self.assertIn('core.neo_dfttest.DFTTest(src, mode="debug", planes=[0], sigma=8.0, tbsize=1)', script)

    def test_renders_avs_filter_call_with_arrays_as_strings(self):
        source = {
            "type": "blank",
            "format": "Y8",
            "width": 64,
            "height": 48,
            "length": 9,
            "color": 64,
        }
        params = {"tbsize": 1, "sigma": 8.0, "planes": [0]}

        script = baseline_runner.render_avs_script("/tmp/plugin.so", source, params)

        self.assertIn('LoadPlugin("/tmp/plugin.so")', script)
        self.assertIn('src = BlankClip(width=64, height=48, length=9, pixel_type="Y8", color=64)', script)
        self.assertIn('return neo_dfttest(src, planes="0", sigma=8.0, tbsize=1)', script)

    def test_renders_avs_yuv_source_color_without_quoting(self):
        source = {
            "type": "blank",
            "format": "Y8",
            "width": 64,
            "height": 48,
            "length": 9,
            "color": 64,
            "avs_color_yuv": "$404040",
        }

        script = baseline_runner.render_avs_script("/tmp/plugin.so", source, {})

        self.assertIn('src = BlankClip(width=64, height=48, length=9, pixel_type="Y8", color_yuv=$404040)', script)

class VapourSynthSourceSelectionTests(unittest.TestCase):
    class SourceNamespace:
        def __init__(self, marker):
            self.marker = marker

        def Source(self, source):
            return ("ffms2", self.marker, source)

        def LibavSMASHSource(self, source):
            return ("lsmas", self.marker, source)

    def test_selects_ffms2_before_lsmas_when_both_are_available(self):
        core = type("Core", (), {})()
        core.ffms2 = self.SourceNamespace("a")
        core.lsmas = self.SourceNamespace("b")
        source = {"type": "ffms2", "resolved_path": "/tmp/source.mp4"}

        clip = baseline_runner._vs_source_clip(core, object(), source)

        self.assertEqual(clip, ("ffms2", "a", "/tmp/source.mp4"))

    def test_falls_back_to_lsmas_when_ffms2_is_not_available(self):
        core = type("Core", (), {})()
        core.lsmas = self.SourceNamespace("b")
        source = {"type": "ffms2", "resolved_path": "/tmp/source.mp4"}

        clip = baseline_runner._vs_source_clip(core, object(), source)

        self.assertEqual(clip, ("lsmas", "b", "/tmp/source.mp4"))

    def test_raises_clear_error_when_no_media_source_filter_is_available(self):
        source = {"type": "ffms2", "resolved_path": "/tmp/source.mp4"}

        with self.assertRaisesRegex(RuntimeError, "no supported VapourSynth source filter found"):
            baseline_runner._vs_source_clip(type("Core", (), {})(), object(), source)


if __name__ == "__main__":
    unittest.main()
