package util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.NoSuchElementException;

import reconstructionAlgorithm.PhylogenyTree;

public class Oracle {

	public List<Quartet> quartetList;
	public HashMap<String, Integer> taxonSet;
	
	public Oracle() {
		quartetList = new ArrayList<Quartet>();
		taxonSet = new HashMap<String, Integer>();
	}
	
	/*
	 * add the given quartet to the Oracle
	 * update fields quartetList and taxonSet
	 */
	public void add(Quartet q) {
		String name[] = new String[4];
		name[0] = q.a.name;
		name[1] = q.b.name;
		name[2] = q.c.name;
		name[3] = q.d.name;
		//update taxonSet by adding keys that do not exist in taxonSet
		for (int i = 0; i < 4; i++) {
			if (taxonSet.containsKey(name[i])) continue;
			else {
				System.out.println("PUT");
				taxonSet.put(name[i], 1);
			}
		}
		quartetList.add(q);
	}
	
	/*
	 * return the quartet topology of a, b, c, d
	 * this function will return the quartet that contains names specified
	 * by the inputs
	 * throw exception if such quartet is not found
	 */
	public Quartet findQuartet(String a, String b, String c, String d) {
		HashMap<String, String> map = new HashMap<String, String>();
		
		for (Quartet i: quartetList) {
			map.put(i.a.name, i.a.name);
			map.put(i.b.name, i.b.name);
			map.put(i.c.name, i.c.name);
			map.put(i.d.name, i.d.name);
			
			if (map.containsKey(a) && map.containsKey(b) &&
				map.containsKey(c) && map.containsKey(d)) {
				return i;
			}
			else {
				map.clear();
				continue;
			}
		}
		throw new NoSuchElementException();
	}
	
	/*
	 * construct the oracle object based on a list of taxons we are going
	 * too study and a systematic PhylogenyTree assumed to be true
	 */
	public void constructOracle (String[] taxons, PhylogenyTree t) {
		int n = taxons.length;
		List<NameSet> list = new ArrayList<NameSet>();
		
		for (int i = 0; i < taxons.length - 3; i++) {
			for (int j = i + 1; j < taxons.length - 2; j++) {
				for (int k = j + 1; k < taxons.length - 1; k++) {
					for (int l = k + 1; l < taxons.length; l++) {
						NameSet s = new NameSet(taxons[i], taxons[j]
								, taxons[k], taxons[l]);
						list.add(s);
					}
				}
			}
		}
		
		//for each NameSet in list, create correct quartet topology and 
		//add the quartet to the Oracle
		for (NameSet q: list) {
			int distance1 = t.findDistance(q.names[0], q.names[1]);
			int distance2 = t.findDistance(q.names[0], q.names[2]);
			int distance3 = t.findDistance(q.names[0], q.names[3]);
			int distance4 = t.findDistance(q.names[1], q.names[2]);
			int distance5 = t.findDistance(q.names[1], q.names[3]);
			int distance6 = t.findDistance(q.names[2], q.names[3]);
			
			int[] distance = {distance1, distance2, distance3, distance4,
					         distance5, distance6};
			System.out.println(q.names[0]+q.names[1]+q.names[2]+q.names[3]);
			System.out.println(distance1 + " " + distance2 + " " + distance3 + " " + distance4 + " " + distance5 + " " + distance6);
			int min = minIndex(distance);
			System.out.println("minDistance: " + min);
			Quartet newquartet = null;
			if (min == 0) {
				newquartet = new Quartet(q.names[0], q.names[1], q.names[2], q.names[3]);
				this.add(newquartet);
			}
			else if (min == 1) {
				newquartet = new Quartet(q.names[0], q.names[2], q.names[1], q.names[3]);
				this.add(newquartet);
			}
			else if (min == 2) {
				newquartet = new Quartet(q.names[0], q.names[3], q.names[1], q.names[2]);
				this.add(newquartet);
			}
			else if (min == 3) {
				newquartet = new Quartet(q.names[0], q.names[3], q.names[1], q.names[2]);
				this.add(newquartet);
			}
			else if (min == 4) {
				newquartet = new Quartet(q.names[0], q.names[2], q.names[1], q.names[3]);
				this.add(newquartet);
			}
			else {
				newquartet = new Quartet(q.names[0], q.names[1], q.names[2], q.names[3]);
				this.add(newquartet);
			}
			System.out.println(newquartet.a.name + " " + newquartet.b.name + " " + 
					newquartet.c.name + " " + newquartet.d.name);
			/*
			int distance1 = t.findDistance(q.names[0], q.names[1]);
			int distance2 = t.findDistance(q.names[0], q.names[2]);
			int distance3 = t.findDistance(q.names[0], q.names[3]);
			int[] distance = {distance1, distance2, distance3};
			
			if (distance1 == distance2 && distance2 == distance3) {
				
			}
			System.out.println(q.names[0]+q.names[1]+q.names[2]+q.names[3]);
			System.out.println(distance1 + " " + distance2 + " " + distance3);
			int min = minIndex(distance);
			System.out.println("minDistance: " + min);
			Quartet newquartet = null;
			if (min == 0) {
				newquartet = new Quartet(q.names[0], q.names[1], q.names[2], q.names[3]);
				this.quartetList.add(newquartet);
			}
			else if (min == 1) {
				newquartet = new Quartet(q.names[0], q.names[2], q.names[1], q.names[3]);
				this.quartetList.add(newquartet);
			}
			else {
				newquartet = new Quartet(q.names[0], q.names[3], q.names[1], q.names[2]);
				this.quartetList.add(newquartet);
			}
			System.out.println(newquartet.a.name + " " + newquartet.b.name + " " + 
					newquartet.c.name + " " + newquartet.d.name);
			*/
		}
	}
	
	/*
	 * return the index of number that contains minimum value across the 
	 * array
	 */
	public static int minIndex(int[] arr) {
		int res = 0;
		for (int i = 1; i < arr.length; i++) {
			if (arr[i] < arr[res]) res = i;
		}
		return res;
	}
	
	/*
	 * a utility object that contains the names of four taxons that 
	 * we are going to study
	 */
	public static class NameSet {
		public String[] names = new String[4];
		
		public NameSet(String a, String b, String c, String d) {
			names[0] = a;
			names[1] = b;
			names[2] = c;
			names[3] = d;
		}
		
		public void printInfo() {
			System.out.println(names[0] + " " + names[1] + " " 
		                       + names[2] + " " + names[3]);
		}
	}
	
	
	public static void main(String[] args) {
		List<String> t = new ArrayList<String>(6);
		t.add("a");
		t.add("b");
		t.add("c");
		t.add("d");
		t.add("e");
		t.add("f");
		PhylogenyTree s = null;
		Oracle Q = new Oracle();
	}
}
