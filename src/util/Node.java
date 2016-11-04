package util;

public class Node {

	
	public String name;
	public Node[] neighbor;
	public Position2D posn = new Position2D();
	public int type; //0 for ancestor nodes and 1 for leaves 2 for super node
	public boolean visited = false;
	/*
	 * conditional constructor
	 * given the name of the node and whether it is a leaf (taxon) or not
	 */
	public Node(String myname, boolean isLeaf) {
		if (isLeaf) {
			name = myname;
			neighbor = new Node[1];
			type = 1;
		}
		else {
			name = myname;
			neighbor = new Node[3];
			type = 0;
		}
	}
	
	/*
	 * add the given node to be this Node's neighbors
	 */
	public void addNeighbors(Node a, Node b, Node c) {
		this.neighbor[0] = a;
		this.neighbor[1] = b;
		this.neighbor[2] = c;
	}
	
	/*
	 * let the given node be the direct ancestor of this leaf node
	 */
	public void setAncestor(Node ances) {
		this.neighbor[0] = ances;
	}
	
	/*
	 * get the type of the node
	 * return true if its a leaf
	 * false otherwise
	 */
	public boolean isLeaf() {
		return this.type == 1;
	}

	/*
	 * return true if the node's type equals 2
	 * false otherwise
	 */
	public boolean isSuperTaxon() {
		return this.type == 2;
	}
}
