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


if __name__ == "__main__":
    unittest.main()
