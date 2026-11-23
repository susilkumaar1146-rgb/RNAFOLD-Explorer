"""
Checks on the Streamlit interface files that do not need Streamlit itself:
every page referenced by the navigation and page links exists and compiles,
and no calculation module imports Streamlit.
"""

import os
import py_compile
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGES = ["home", "sequence_analysis", "structure_prediction", "mutation_analysis",
         "visualization", "results_downloads"]


def read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as handle:
        return handle.read()


def test_all_pages_exist_and_compile():
    for page in PAGES:
        path = os.path.join(ROOT, "app_pages", page + ".py")
        assert os.path.exists(path), path
        py_compile.compile(path, doraise=True)
    for path in ("app.py", "ui/state.py", "ui/components.py"):
        py_compile.compile(os.path.join(ROOT, path), doraise=True)


def test_navigation_and_links_point_to_real_pages():
    referenced = set(re.findall(r'"(app_pages/[a-z_]+\.py)"', read("app.py") + read("ui/components.py")))
    assert referenced, "no page paths found"
    for path in referenced:
        assert os.path.exists(os.path.join(ROOT, path)), path
    assert len(re.findall(r"st\.Page\(", read("app.py"))) == len(PAGES)


def test_calculation_modules_do_not_import_streamlit():
    for filename in os.listdir(os.path.join(ROOT, "modules")):
        if filename.endswith(".py"):
            assert "import streamlit" not in read(os.path.join("modules", filename)), filename


def test_widget_keys_unique_per_page():
    for page in PAGES:
        keys = re.findall(r'key="([^"{]+)"', read(os.path.join("app_pages", page + ".py")))
        assert len(keys) == len(set(keys)), (page, keys)
