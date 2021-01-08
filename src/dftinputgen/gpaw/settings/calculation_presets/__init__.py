import os
import json
import pkg_resources


__all__ = ["GPAW_PRESETS"]


GPAW_PRESETS = {}


preset_listdir = pkg_resources.resource_listdir(
    "dftinputgen.gpaw.settings", "calculation_presets"
)
for filename in preset_listdir:
    root, ext = os.path.splitext(filename)
    if not ext == ".json":
        continue
    resource = pkg_resources.resource_filename(
        "dftinputgen.gpaw.settings.calculation_presets", filename
    )
    with open(resource, "r") as fr:
        GPAW_PRESETS[root] = json.load(fr)
