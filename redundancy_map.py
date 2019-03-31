def add_to_map(seq, map, end_size): #required nested map
	## if key does not exist, add it
	new_key = seq[-20:]
	if not map[new_key]:

	#new_value = {"sequence": seq, "length": len(seq), "depth": 0}
	map.add()

if __name__ == "__main__":
	## stores non-redundant integrations in nested dictionary
	nested_dict = lambda: collections.defaultdict(nested_dict)
	redundancy_map = nested_dict()
