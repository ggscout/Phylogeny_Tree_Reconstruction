package util;

public class Quartet {

	//four leaf nodes
	public Node a = new Node("", true);
	public Node b = new Node("", true);
	public Node c = new Node("", true);
	public Node d = new Node("", true);
	//two ancestor nodes
	public Node ancesone = new Node("ances", false);
	public Node ancestwo = new Node("ances", false); 	
	
	/*
	 * construct an arbitrary quartet topology given the four names in 
	 * topology order
	 * Notice that the name MUST BE given in topology order
	 */
	public Quartet(String namea, String nameb, String namec, String named) {
		a.neighbor[0] = ancesone;
		b.neighbor[0] = ancesone;
		ancesone.neighbor[0] = a;
		ancesone.neighbor[1] = b;
		ancesone.neighbor[2] = ancestwo;
		ancestwo.neighbor[2] = ancesone;
		ancestwo.neighbor[0] = c;
		ancestwo.neighbor[1] = d;
		c.neighbor[0] = ancestwo;
		d.neighbor[0] = ancestwo;
		a.name = namea;
		b.name = nameb;
		c.name = namec;
		d.name = named;
	}
	
	/*
	 * return the name of the node which stands at the same side
	 * with the node named as x
	 */
	public String findConnection(String x) {
		if (a.name.equals(x)) {
			return b.name;
		}
		else if (b.name.equals(x)) {
			return a.name;
		}
		else if (c.name.equals(x)) {
			return d.name;
		}
		else {
			return c.name;
		}
	}
}
