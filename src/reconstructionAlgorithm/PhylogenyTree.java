package reconstructionAlgorithm;
import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Stack;
import util.Node;
import util.Position2D;
import util.Oracle;
import util.Quartet;
import util.StdDraw;

public class PhylogenyTree {

	//represents all leaves from the tree using their names
	public ArrayList<String> taxons = new ArrayList<String>();
	//represents all the nodes in the PhylogenyTree
	public ArrayList<Node> nodes = new ArrayList<Node>();
	
	
	//auxiliary fields for Algorithm usage
	public List<Node> leaves = new ArrayList<Node>();
	//auxiliary fields for Algorithm usage
	public List<Node> separators = new ArrayList<Node>();
	//auxiliary fields for Algorithm usage
	public Node currSepara = null;
	//auxiliary fields for Algorithm usage
	public HashMap<Node, Integer> visited = new HashMap<Node, Integer>();
	
	
	/*
	 * default constructor
	 */
	public PhylogenyTree() {
		taxons = new ArrayList<String>();
		nodes = new ArrayList<Node>();
	}
	
	/*
	 * constructor
	 */
	public PhylogenyTree(String[] t, Node[] n) {
		for (String i: t) {
			this.taxons.add(i);
		}
		for (Node j: n) {
			this.nodes.add(j);
		}
	}
	
	/*
	 * given leaves a, b, c, d, x and oracle Q
	 * return the "true" quartet topology over a, b, c, d
	 */
	public static Quartet baseQuery(String a, String b, String c, String d, String x, Oracle Q) {
		Quartet q1 = Q.findQuartet(a, b, c, x);
		String y1, y2, y3;
		y3 = q1.findConnection(x);
		if (y3.equals(c)) {
			y1 = a;
			y2 = b;
		}
		else if (y3.equals(b)) {
			y1 = a;
			y2 = c;
		}
		else {
			y1 = b;
			y2 = c;
		}
		
		Quartet q2 = Q.findQuartet(y1, y2, x, d);
		
		Quartet q3;
		if (q2.findConnection(x).equals(d)) {
			q3 = new Quartet(y1, y2, y3, d);
			return q3;
		}
		else if (q2.findConnection(x).equals(y1)) {
			q3 = new Quartet(y1, y3, y2, d);
			return q3;
		}
		else {
			q3 = new Quartet(y1, d, y2, y3);
			return q3;
		}
	}
	
	/*
	 * Given quartet set Oracle Q
	 * return the most accurate topology on a, b, c, d by running BaseQuery
	 * on every node x in Q - {a, b, c, d}
	 */
	public static Quartet highProbQuery(String a, String b, String c, String d, Oracle Q) {
		int count[] = new int[3]; //ab|cd ac|bd ad|bc
		for (int i = 0; i < 3; i++) {
			count[i] = 0;
		}
		String x = "";
		for (String key: Q.taxonSet.keySet()) {
			if (key.equals(a) || key.equals(b)
				|| key.equals(c) || key.equals(d)) {
				continue;
			}
			else {
				x = key;
				Quartet q = baseQuery(a, b, c, d, x, Q);
				if (q.findConnection(a).equals(b)) {
					count[0]++;
				}
				else if (q.findConnection(a).equals(c)) {
					count[1]++;
				}
				else {
					count[2]++;
				}
			}
		}
		
		int maxIndex = findMaxIndex(count);
		if (maxIndex == 0) {
			Quartet result = new Quartet(a, b, c, d);
			return result;
		}
		else if (maxIndex == 1) {
			Quartet result = new Quartet(a, c, b, d);
			return result;
		}
		else {
			Quartet result = new Quartet(a, d, b, c);
			return result;
		}	
	}
	

	/*
	 * PhylogenyTree Reconstruction function
	 * update the "nodes" field of this PhylogenyTree object
	 * the "taxons" field of this object is already given because we know
	 * what biological species (taxa) we are studying at
	 */
	public void highProbReconstruct(Oracle Q) {
		//start with any four taxa
		String a, b, c, d;
		a = taxons.get(0);
		b = taxons.get(1);
		c = taxons.get(2);
		d = taxons.get(3);
		Quartet q = highProbQuery(a, b, c, d, Q);
		nodes.add(q.a);
		nodes.add(q.b);
		nodes.add(q.c);
		nodes.add(q.d);
		nodes.add(q.ancesone);
		nodes.add(q.ancestwo);
		for (int i = 4; i < taxons.size(); i++) {
			Node s = this.findSeparator();
			insert(taxons.get(i), Q, s, null);
		}
		
		
		//updating the PhylogenyTree 
		//pick a random ancestor node first
		Node r = null;
		for (int i = 0; i < this.nodes.size(); i++) {
			if (this.nodes.get(i).isLeaf())  {
				this.visited.put(this.nodes.get(i), 1);
				continue;
			}
			else {
				this.visited.put(this.nodes.get(i), 1);
				r = this.nodes.get(i);
			}
		}
		
		PhylogenyTree updated = new PhylogenyTree();
		updated.DFS(r);
		for (Node i: updated.nodes) {
			if (!this.visited.containsKey(i)) {
				this.nodes.add(i);
			}
		}
	}

	/*
	 * insert the node specified by its name s to the PhylogenyTree
	 */
	public void insert(String s, Oracle Q, Node separator, Node prevSepara) {
		//random array contains random leaf nodes from three branches of 
		//the separator
		Node[] random = new Node[3];
		for (int i = 0; i < 3; i++) {
			if (separator.neighbor[i].isLeaf() || separator.neighbor[i].isSuperTaxon()) {
				random[i] = separator.neighbor[i];
			}
			else {
				random[i] = findRandomLeaf(separator.neighbor[i], separator);
			}
		}
		
		String[] names = new String[3];
		for (int i = 0; i < 3; i++) {
			if (random[i].isLeaf()) {
				names[i] = random[i].name;
			}
			else {
				int nameLength = random[i].name.length();
				int ranIndex = (int) (Math.random() * nameLength);
				names[i] = random[i].name.substring(ranIndex, ranIndex + 1);
			}			
		}
		
		Quartet q = highProbQuery(names[0], names[1], names[2], s, Q);
		Node connectingNode = random[0];
		int index = 0;
		for (int i = 0; i < 3; i++) {
			if (q.findConnection(s).equals(names[i])) {
				connectingNode = random[i];
				index = i;
				break;
			}
		}
		
		if (connectingNode.isLeaf() && connectingNode.neighbor[0] == separator) {
			Node newNode = new Node(s, true);
			Node newAnces = new Node("ances", false);
			
			
			newAnces.neighbor[0] = connectingNode;
			connectingNode.neighbor[0] = newAnces;
			
			newAnces.neighbor[1] = newNode;
			newNode.neighbor[0] = newAnces;
			
			newAnces.neighbor[2] = separator;
			separator.neighbor[index] = newAnces;
			
			this.nodes.add(newNode);
			this.nodes.add(newAnces);
			this.separators.clear();
			if (prevSepara != null) {
				Node c = null;
				for (int i = 0; i < 3; i++) {
					if (prevSepara.neighbor[i].isLeaf()) continue;
					else {
						c = prevSepara.neighbor[i];
						for (int j = 0; j < 3; j++) {
							if (c.neighbor[j].isSuperTaxon()) {
								c.neighbor[j].neighbor[0] = null;
								c.neighbor[j] = prevSepara;
								break;
							}
						}
					}
				}
			}
			
			return;
		}
		else if (connectingNode.isSuperTaxon() && connectingNode.neighbor[0] == separator) {
			Node newNode = new Node(s, true);
			Node newAnces = new Node("ances", false);
			
			connectingNode.neighbor[0] = null;
			connectingNode = null;
			
			newAnces.neighbor[0] = prevSepara;
			for (int i = 0; i < 3; i++) {
				if (prevSepara.neighbor[i] == separator) {
					prevSepara.neighbor[i] = newAnces;
				}
			}
			
			newAnces.neighbor[1] = newNode;
			newNode.neighbor[0] = newAnces;
			
			newAnces.neighbor[2] = separator;
			separator.neighbor[index] = newAnces;
			
			this.nodes.add(newNode);
			this.nodes.add(newAnces);
			this.separators.clear();
			return;
		}
		else {
			Node superTaxon = new Node("", true);
			superTaxon.type = 2;
			findAllLeaves(separator, separator.neighbor[index]); 
			for (Node i: leaves) {
				superTaxon.name = superTaxon.name + i.name;
			}
			leaves.clear();
			
			if (prevSepara != null) {
				Node c = null;
				for (int i = 0; i < 3; i++) {
					if (prevSepara.neighbor[i].isLeaf()) continue;
					else {
						c = prevSepara.neighbor[i];
						for (int j = 0; j < 3; j++) {
							if (c.neighbor[j].isSuperTaxon()) {
								c.neighbor[j].neighbor[0] = null;
								c.neighbor[j] = prevSepara;
								break;
							}
						}
					}
				}
			}
			
			Node prev = separator;

			for (int i = 0; i < 3; i++) {
				if (separator.neighbor[index].neighbor[i] == separator) {
					separator.neighbor[index].neighbor[i] = superTaxon;
					superTaxon.neighbor[0] = separator.neighbor[index];
					break;
				}
			}
			
			//create a new PhylogenyTree with this supertaxon
			PhylogenyTree newTree = new PhylogenyTree();
			newTree.DFS(superTaxon.neighbor[0]);
			newTree.visited.clear();
			Node newSeparator = newTree.findSeparator();
			newTree.insert(s, Q, newSeparator, prev);
			return;
		}
	}
	
	/*
	 * return an arbitrary leaf node in the subtree induced by the input 
	 * node; return the input node itself if the subtree is a leaf
	 * else return a random leaf in the subtree starting with node as
	 * the root and excluding the separator
	 */
	public Node findRandomLeaf(Node node, Node separator) {
		if (node.isLeaf()) return node;
		else {
			findAllLeaves(node, separator);
			int randomIndex = (int) (Math.random() * leaves.size());
			Node res = leaves.get(randomIndex);
			leaves.clear();
			return res;
		}
		
	}
	
	/*
	 * find all the leaves in the subtree starting with node as the root
	 * and excluding the separator
	 * adding all the leaves to instance variable list leaves 
	 */
	public void findAllLeaves(Node node, Node separator) {
		for (int i = 0; i < 3; i++) {
			if (node.neighbor[i] == separator) continue;
			else if (node.neighbor[i].isLeaf() || node.neighbor[i].isSuperTaxon()) {
				leaves.add(node.neighbor[i]);
			}
			else {
				findAllLeaves(node.neighbor[i], node);
			}
		}
	}
	
	/*
	 * return a random valid separator node from the PhylogenyTree
	 * take visited HashMap as input, which represents the separator 
	 * already visited during this round of insertion
	 */
	public Node findSeparator() {
		this.initializeSeparators();
		int index = (int) (Math.random() * this.separators.size());
		return this.separators.get(index);
	}
	
	/*
	 * initialize the empty separators list with all the internal nodes
	 * delete those nodes that are no eligible to become a separator
	 */
	public void initializeSeparators() {
		if (!this.separators.isEmpty()) return;
		
		//add all internal nodes
		for (int i = 0; i < this.nodes.size(); i++) {
			if (nodes.get(i).isLeaf() || nodes.get(i).isSuperTaxon()) continue;
			else {
				this.separators.add(this.nodes.get(i));
			}
		}
		
		//throw exceptions
		if (this.separators.isEmpty()) throw new NoSuchElementException();
		//remove invalid internal nodes
		for (int i = 0; i < this.separators.size(); i++) {
			Node n = this.separators.get(i);
			if (!isValidSeparator(n)) {
				this.separators.remove(i);
				i--;
			}
		}
	}
	/*
	 * utility function
	 * return true if the node x in the PhylogenyTree is a valid separator
	 * each subtree generated by picking x as a separator must have at most
	 * half of the leaves of the PhylogenyTree
	 * return false if it is not valid
	 */
	public static boolean isValidSeparator(Node x) {
		int[] count = new int[3];
		for (int i = 0; i < 3; i++) {
			count[i] = countLeaves(x.neighbor[i], x);
		}
		int maxIndex = findMaxIndex(count);
		int totalLeaves = count[0] + count[1] + count[2];
		int half = totalLeaves / 2;
		if (count[maxIndex] <= half) return true;
		else return false;
	}
	
	/*
	 * utility function
	 * return the number of leaves of the subtree containing x
	 * the subtree is created by picking the separator
	 */
	public static int countLeaves(Node x, Node separator) {
		if (x.isLeaf() || x.isSuperTaxon()) return 1;
		else {
			int count = 0;
			for (int i = 0; i < 3; i++) {
				if (x.neighbor[i] == separator) continue;
				else {
					if (x.neighbor[i].isLeaf()) count++;
					else {
						count += countLeaves(x.neighbor[i], x);
					}
				}
			}
			return count;
		}
	}
	
	/*
	 * starting with a ancestor node
	 * add all leaf & ancestor nodes from this connected component to nodes
	 * list
	 */
	public void DFS(Node root) {
		this.nodes.add(root);
		if (root.isLeaf() || root.isSuperTaxon()) return;
		else {
			for (int i = 0; i < 3; i++) {
				if (root.neighbor[i] == null) continue;
				else if (!visited.containsKey(root.neighbor[i])) {
					visited.put(root.neighbor[i], 1);
					DFS(root.neighbor[i]);
				}
				else continue;
			}
		}
		
	}

	/*
	 * utility function
	 * return the index for the maximum value in an array
	 */
	public static int findMaxIndex(int[] a) {
		int maxIndex = 0;
		for (int i = 0; i < a.length; i++) {
			if (a[i] > a[maxIndex]) {
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	/*
	 * find the distance between the given root node of the tree and its
	 * farthest leaf; distance is denoted as the number of edges between
	 * these two nodes
	 */
	public static int findMaxDepth(Node root) {
		root.visited = true;
		if (root.isLeaf()) 
		{
			return 0;
		}
		else {
			int d0, d1, d2;
			d0 = 0;
			d1 = 0;
			d2 = 0;
			if (!root.neighbor[0].visited) {
				root.neighbor[0].visited = true;
				d0 = findMaxDepth(root.neighbor[0]);
			}
			if (!root.neighbor[1].visited) {
				root.neighbor[1].visited = true;
				d1 = findMaxDepth(root.neighbor[1]);
			}
			
			if (!root.neighbor[2].visited) {
				root.neighbor[2].visited = true;
				d2 = findMaxDepth(root.neighbor[2]);
			}
			return maxofThree(d0, d1, d2) + 1;
		}
	}
	
	/*
	 * return the max number of the given three numbers
	 */
	public static int maxofThree(int a, int b, int c) {
		int[] arr = {a, b, c};
		int max = a;
		for (int i = 0; i < 3; i++) {
			if (arr[i] > a) max = arr[i];
		}
		return max;
	}
	/*
	 * print out relevant information of the PhylogenyTree
	 * including the name of its nodes and their neighbors
	 */
	public static void printTreeInfo(PhylogenyTree x) {
		System.out.println("There are a total of " + x.nodes.size() + " animal/plant species");
		List<Node> ances = new ArrayList<Node>();
		List<Node> leaves = new ArrayList<Node>();
		for (int i = 0; i < x.nodes.size(); i++) {
			if (!x.nodes.get(i).isLeaf()) {
				ances.add(x.nodes.get(i));
			}
			else {
				leaves.add(x.nodes.get(i));
			}
		}
		for (int i = 0; i < ances.size(); i++) {
			String s = Integer.toString(i+1);
			ances.get(i).name = s;
		}
		
		System.out.println("THESE ARE THE ANCESTOR ANIMAL/PLANT SPECIES:");
		for (int i = 0; i < ances.size(); i++) {
			System.out.println("Here is Ancestor No. " + ances.get(i).name);
			System.out.println("Its neighbors are: ");
			for (int j = 0; j < 3; j++) {
				System.out.println(ances.get(i).neighbor[j].name);
			}
		}
		
		System.out.println("THESE ARE THE TAXONS:");
		for (int i = 0; i < leaves.size(); i++) {
			System.out.println("Here is taxon: " + leaves.get(i).name);
			for (int j = 0; j < leaves.get(i).neighbor.length; j++) {
				System.out.println("Its direct ancestor is: " + leaves.get(i).neighbor[j].name);
			}
		}
	}
	
	
	/*
	 * draw the PhylogenyTree using StdDraw GUI provided by Princeton
	 * Professor Robert Sedgewick
	 */
	/*
	public void drawTree(Node root) {
		int depth = findMaxDepth(root);
		this.setNodeUnvisited();
		root.posn = new Position2D(50.0, 50.0);
		root.visited = true;
		root.neighbor[0].posn = new Position2D(50.0, 70.0);
		double rad = Math.toRadians(30);
		root.neighbor[1].posn = new Position2D(50.0 - 20.0 * Math.cos(rad), 40.0);
		root.neighbor[2].posn = new Position2D(50.0 + 20.0 * Math.cos(rad), 40.0);
		drawTreeHelper(root.neighbor[0]);
		drawTreeHelper(root.neighbor[1]);
		drawTreeHelper(root.neighbor[2]);
	}
	
	
	public void drawTreeHelper(Node n) {
		n.visited = true;
		if (n.isLeaf()) return;
		else {
			Position2D center = n.posn;
			Position2D visitedneighbor = null;
			for (int i = 0; i < 3; i++) {
				if (n.neighbor[i].visited) {
					visitedneighbor = n.neighbor[i].posn;
					break;
				}
			}
			
		}
	}
	*/
	
	
	/*
	 * set the visited field of all nodes to be false
	 */
	public void setNodeUnvisited() {
		for (Node i: this.nodes) {
			i.visited = false;
		}
	}
	
	/*
	 * find the distance (number of edges on the path) between two leaves
	 * named x and y
	 */
	public int findDistance(String x, String y) {
		Node a = null;
		this.setNodeUnvisited();
		for (int i = 0; i < this.nodes.size(); i++) {
			if (this.nodes.get(i).name.equals(x)) {
				a = this.nodes.get(i);
				break;
			}
		}
		int count = 1;
		a.visited = true;
		return findDistanceHelper(a, y, count);
	}
	
	/*
	 * use breath first search to find the distance between a leaf node n 
	 * and another leaf node named y
	 */
	public int findDistanceHelper(Node n, String y, int count) {
		Queue<Node> q = new LinkedList<Node>();
		Queue<Node> level = new LinkedList<Node>();
	    q.add(n.neighbor[0]);
		level.add(n.neighbor[0]);
		while (!q.isEmpty()) {
			int queueSize = q.size();
			for (int i = 0; i < queueSize; i++) {
				Node d = q.poll();
				d.visited = true;
				if (d.name.equals(y)) {
					return count;
				}
			}		
			for (Node j: level) {
				if (j.isLeaf()) continue;
				for (int k = 0; k < 3; k++) {
					if (!j.neighbor[k].visited) {
						q.add(j.neighbor[k]);
					}
				}
			}
			level.clear();
			for (Node m: q) {
				level.add(m);
			}
			count++;
		}
		return count;
	}
	
	/*
	 * set the internal nodes' names
	 */
	public void setAncesName() {
		List<Node> l = new ArrayList<Node>();
		for (int i = 0; i < this.nodes.size(); i++) {
			if (this.nodes.get(i).isLeaf()) continue;
			else {
				l.add(this.nodes.get(i));
			}
		}
		for (int j = 0; j < l.size(); j++) {
			l.get(j).name = Integer.toString(j + 1);
		}
	}
	
	public static void main(String[] args) {
		//Test Case 1
		/*
		PhylogenyTree x = new PhylogenyTree();
		Quartet q1 = new Quartet("a", "b", "c", "d");
		Quartet q2 = new Quartet("a", "b", "c", "e");
		Quartet q3 = new Quartet("a", "b", "c", "f");
		Quartet q4 = new Quartet("a", "b", "d", "e");
		Quartet q5 = new Quartet("a", "b", "d", "f");
		Quartet q6 = new Quartet("a", "b", "e", "f");
		Quartet q7 = new Quartet("a", "c", "d", "e");
		Quartet q8 = new Quartet("a", "c", "d", "f");
		Quartet q9 = new Quartet("a", "c", "e", "f");
		Quartet q10 = new Quartet("a", "d", "e", "f");
		Quartet q11 = new Quartet("b", "c", "d", "e");
		Quartet q12 = new Quartet("b", "c", "d", "f");
		Quartet q13 = new Quartet("b", "c", "e", "f");
		Quartet q14 = new Quartet("b", "d", "e", "f");
		Quartet q15 = new Quartet("c", "d", "e", "f");
		Oracle Q = new Oracle();
		Q.add(q1);
		Q.add(q2);
		Q.add(q3);
		Q.add(q4);
		Q.add(q5);
		Q.add(q6);
		Q.add(q7);
		Q.add(q8);
		Q.add(q9);
		Q.add(q10);
		Q.add(q11);
		Q.add(q12);
		Q.add(q13);
		Q.add(q14);
		Q.add(q15);
		x.taxons.add("a");
		x.taxons.add("b");
		x.taxons.add("c");
		x.taxons.add("d");
		x.taxons.add("e");
		x.taxons.add("f");
		x.highProbReconstruct(Q);
		x.setAncesName();
		Node root = null;
		for (int i = 0; i < 10; i++) {
			if (!x.nodes.get(i).isLeaf()) {
				root = x.nodes.get(i);
				break;
			}
		}
		System.out.println(x.findDistance("a", "f"));
		*/
		
		/*
		//Test Case 2
		PhylogenyTree y = new PhylogenyTree();
		Quartet n1 = new Quartet("a", "b", "c", "d");
		Quartet n2 = new Quartet("a", "b", "c", "e");
		Quartet n3 = new Quartet("a", "b", "c", "f");
		Quartet n4 = new Quartet("a", "b", "d", "e");
		Quartet n5 = new Quartet("a", "b", "d", "f");
		Quartet n6 = new Quartet("a", "b", "e", "f");
		Quartet n7 = new Quartet("a", "e", "c", "d");
		Quartet n8 = new Quartet("a", "f", "c", "d");
		Quartet n9 = new Quartet("a", "c", "e", "f");
		Quartet n10 = new Quartet("a", "d", "e", "f");
		Quartet n11 = new Quartet("b", "e", "c", "d");
		Quartet n12 = new Quartet("b", "f", "c", "d");
		Quartet n13 = new Quartet("b", "c", "e", "f");
		Quartet n14 = new Quartet("b", "d", "e", "f");
		Quartet n15 = new Quartet("c", "d", "e", "f");
		Oracle Q1 = new Oracle();
		Q1.add(n1);
		Q1.add(n2);
		Q1.add(n3);
		Q1.add(n4);
		Q1.add(n5);
		Q1.add(n6);
		Q1.add(n7);
		Q1.add(n8);
		Q1.add(n9);
		Q1.add(n10);
		Q1.add(n11);
		Q1.add(n12);
		Q1.add(n13);
		Q1.add(n14);
		Q1.add(n15);
		y.taxons.add("a");
		y.taxons.add("b");
		y.taxons.add("c");
		y.taxons.add("d");
		y.taxons.add("e");
		y.taxons.add("f");
		y.highProbReconstruct(Q1);
		Node root = null;
		for (int j = 0; j < y.nodes.size(); j++) {
			if (!y.nodes.get(j).isLeaf()) {
				root = y.nodes.get(j);
				break;
			}
		}
		printTreeInfo(y);
		*/
		
		/*
		//Test Case 3
		Node a = new Node("a", true);
		Node b = new Node("b", true);
		Node c = new Node("c", true);
		Node d = new Node("d", true);
		Node e = new Node("e", true);
		Node f = new Node("f", true);
		Node one = new Node("1", false);
		Node two = new Node("2", false);
		Node three = new Node("3", false);
		Node four = new Node("4", false);
		a.setAncestor(one);
		b.setAncestor(one);
		c.setAncestor(two);
		d.setAncestor(two);
		e.setAncestor(three);
		f.setAncestor(three);
		one.addNeighbors(a, b, four);
		two.addNeighbors(c, d, four);
		three.addNeighbors(e, f, four);
		four.addNeighbors(one, two, three);
		String[] t = {"a", "b", "c", "d", "e", "f"};
		Node[] n = {a, b, c, d, e, f, one, two, three, four};
		PhylogenyTree z = new PhylogenyTree(t, n);
		Oracle Q = new Oracle();
		Q.constructOracle(t, z);
		
		PhylogenyTree l = new PhylogenyTree();
		l.taxons.add("a");
		l.taxons.add("b");
		l.taxons.add("c");
		l.taxons.add("d");
		l.taxons.add("e");
		l.taxons.add("f");
		for (int m = 0; m < Q.quartetList.size(); m++) {
			Quartet q = Q.quartetList.get(m);
			System.out.println(q.a.name+q.b.name+q.c.name+q.d.name);
		}
		for (String km: Q.taxonSet.keySet()) {
			System.out.println(km);
		}
		l.highProbReconstruct(Q);
		printTreeInfo(l);
		*/
		
		/*
		//Test case 3
		Node a = new Node("a", true);
		Node b = new Node("b", true);
		Node c = new Node("c", true);
		Node d = new Node("d", true);
		Node e = new Node("e", true);
		Node f = new Node("f", true);
		Node one = new Node("1", false);
		Node two = new Node("2", false);
		Node three = new Node("3", false);
		Node four = new Node("4", false);
		a.setAncestor(one);
		b.setAncestor(one);
		c.setAncestor(four);
		d.setAncestor(four);
		e.setAncestor(three);
		f.setAncestor(two);
		one.addNeighbors(a, b, two);
		two.addNeighbors(one, three, f);
		three.addNeighbors(two, four, e);
		four.addNeighbors(c, d, three);
		String[] t = {"a", "b", "c", "d", "e", "f"};
		Node[] n = {a, b, c, d, e, f, one, two, three, four};
		PhylogenyTree z = new PhylogenyTree(t, n);
		Oracle Q = new Oracle();
		Q.constructOracle(t, z);
		PhylogenyTree l = new PhylogenyTree();
		l.taxons.add("a");
		l.taxons.add("b");
		l.taxons.add("c");
		l.taxons.add("d");
		l.taxons.add("e");
		l.taxons.add("f");
		l.highProbReconstruct(Q);
		printTreeInfo(l);
		*/
		
	}
}
