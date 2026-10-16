import json
from path_homology import Digraph, compute_h1, verify_certificate

result = compute_h1(Digraph(3, [(0, 1), (1, 2), (2, 0)]), characteristic=0)
assert result.betti == 1
record = result.certificate()
assert verify_certificate(record)
print(json.dumps(record, indent=2))
