package util;

public class Position2D {
	public double x;
	public double y;
	
	public Position2D() {
		x = 0;
		y = 0;
	}
	public Position2D(double a, double b) {
		x = a;
		y = b;
	}
	
	public double getX() {
		return this.x;
	}
	public double getY() {
		return this.y;
	}
}
