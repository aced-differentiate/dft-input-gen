import os
import glob
import json


__all__ = ['GPAW_BASE_RECIPES']


recipe_jsons = glob.glob(os.path.join(os.path.dirname(__file__), '*.json'))

GPAW_BASE_RECIPES = {}
for recipe_json in recipe_jsons:
    recipe_name = os.path.splitext(os.path.basename(recipe_json))[0]
    with open(recipe_json, 'r') as fr:
        GPAW_BASE_RECIPES[recipe_name] = json.load(fr)
