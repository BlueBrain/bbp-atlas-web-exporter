import copy

def flattenTree(tree, id_prop_name="id", children_prop_name="children"):
    """
    Transforms a nested tree into a flat dictionnary where the key is the id
    (with the name given by id_prop_name) and where the children are not objects but a
    list of ids (under the property name given by children_prop_name).
    This creates the properties: "_parent", "_descendants" and "_ascendants"
    """

    PARENT_PROP_NAME = "_parent"
    DESCENDANTS_PROP_NAME = "_descendants"
    ASCENDANTS_PROP_NAME = "_ascendants"

    # let's not damage the original tree by making a deep copy
    tree_copy = copy.deepcopy(tree)

    root_node = tree_copy
    root_node[PARENT_PROP_NAME] = None
    flat_tree = {}
    node_to_explore = [root_node]

    while len(node_to_explore):
        node = node_to_explore.pop()
        node_id = node[id_prop_name]
        flat_tree[node_id] = node

        children_ids = []

        # adding the children to the list of nodes to explore later
        if children_prop_name in node:
            for child in node[children_prop_name]:
                children_ids.append(child[id_prop_name])
                node_to_explore.append(child)
            # is it faster to do:
            # node_to_explore = node_to_explore + node[children_prop_name]

        # replacing the list of node object with a list of node ids
        node[children_prop_name] = children_ids


    # Create a list of descendants that contains all the children and children
    # of children up to the leaf nodes
    for node_id in flat_tree:
        node = flat_tree[node_id]
        node[DESCENDANTS_PROP_NAME] = []
        node_to_explore = [] + node[children_prop_name]

        while len(node_to_explore):
            child_id = node_to_explore.pop()
            node[DESCENDANTS_PROP_NAME].append(child_id)
            child_node = flat_tree[child_id]
            child_node[PARENT_PROP_NAME] = node_id
            node_to_explore = node_to_explore + child_node[children_prop_name]

    # Create a list of ascendants, that contains the parent, grand parent and so on
    # up to the root node.
    # This is convenient for checking if a node is contained by another with more than
    # one level of hierarchy.
    for node_id in flat_tree:
        node = flat_tree[node_id]
        node[ASCENDANTS_PROP_NAME] = []
        node_to_explore = node[PARENT_PROP_NAME]

        while node_to_explore:
            node[ASCENDANTS_PROP_NAME].append(node_to_explore)
            node_to_explore = flat_tree[node_to_explore][PARENT_PROP_NAME]

    return flat_tree
