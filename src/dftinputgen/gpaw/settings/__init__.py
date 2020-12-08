import json
import pkg_resources

__all__ = ['GPAW_TAGS']

tags_file = pkg_resources.resource_filename('dftinpgen.gpaw.settings',
                                            'tags_and_groups.json')

with open(tags_file, 'r') as fr:
    GPAW_TAGS = json.load(fr)
