import json
import pytest
import rdflib

from blue_brain_atlas_web_exporter.json_to_jsonld import hierarchy_json_to_jsonld, TYPE_FOR_FRAMING
from jsonpath_ng import parse
from rdflib import RDF


def test_hierarchy_json_to_jsonld():
    jsoncontent = json.loads(open("./tests/data/test_for_1.json", "r").read())
    jsonld = hierarchy_json_to_jsonld(jsoncontent)
    assert jsonld is not None

    jsonld_graph = rdflib.Graph()
    jsonld_graph.parse(data=json.dumps(jsonld), format="json-ld")

    jsonpath_expr = parse('$..children.[*]')
    matches = jsonpath_expr.find(jsoncontent)

    assert "defines" in jsonld
    assert type(jsonld["defines"]) == list
    assert len(jsonld["defines"]) == len(matches) + 1 # all children + root

    jsonld_hasPart_expr = parse('$..hasPart.[*]')
    jsonld_matches = jsonld_hasPart_expr.find(jsonld)
    assert len(jsonld_matches) == len(matches) + 1 # all children + root

    assert "nsg:tempDefines" not in jsonld
    assert "@context" in jsonld
    assert "@context" in jsonld

    typed_with_frame = jsonld_graph.triples((None, RDF.type, TYPE_FOR_FRAMING))
    assert list(typed_with_frame) is not None
    assert len(list(typed_with_frame)) == 0

