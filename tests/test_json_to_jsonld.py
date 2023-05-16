import os, json
from rdflib import Graph, RDF
from jsonpath_ng import parse

from src.blue_brain_atlas_web_exporter.json_to_jsonld import hierarchy_json_to_jsonld, TYPE_FOR_FRAMING

test_folder = os.environ["TEST_FOLDER"]
input_hierarchy = "test_for_1.json"

def test_hierarchy_json_to_jsonld():
    input_hierarchy_path = os.path.join(test_folder, "data" ,input_hierarchy)
    jsoncontent = json.loads(open(input_hierarchy_path, "r").read())
    try:
        jsonld = hierarchy_json_to_jsonld(jsoncontent)
    except KeyError as key:
        print(f"{key} in {input_hierarchy_path}.\nTime to update the test file?")
        return
    assert jsonld is not None

    jsonld_graph = Graph()
    jsonld_graph.parse(data=json.dumps(jsonld), format="json-ld")

    jsonpath_expr = parse('$..children.[*]')
    matches = jsonpath_expr.find(jsoncontent)

    assert "defines" in jsonld
    assert type(jsonld["defines"]) == list

    jsonld_hasPart_expr = parse('$..hasPart.[*]')
    jsonld_matches = jsonld_hasPart_expr.find(jsonld)
    assert len(jsonld_matches) == len(matches) + 1 # all children + root

    assert "nsg:tempDefines" not in jsonld
    assert "@context" in jsonld
    assert "@context" in jsonld

    typed_with_frame = jsonld_graph.triples((None, RDF.type, TYPE_FOR_FRAMING))
    assert list(typed_with_frame) is not None
    assert len(list(typed_with_frame)) == 0

