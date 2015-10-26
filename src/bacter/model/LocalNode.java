package bacter.model;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class LocalNode {

    public final LocalNode left;
    public final LocalNode right;

    public final double height;

    final int cachedHashCode;

    public LocalNode(LocalNode left, LocalNode right, double height) {
        this.left = left;
        this.right = right;
        this.height = height;

        this.cachedHashCode = computeHashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        if (this.hashCode() != o.hashCode())
            return false;

        LocalNode localNode = (LocalNode) o;

        if (Double.compare(localNode.height, height) != 0) return false;
        if (left != null ? !left.equals(localNode.left) : localNode.left != null)
            return false;
        return !(right != null ? !right.equals(localNode.right) : localNode.right != null);

    }

    public int computeHashCode() {
        int result;
        long temp;
        result = left != null ? left.hashCode() : 0;
        result = 31 * result + (right != null ? right.hashCode() : 0);
        temp = Double.doubleToLongBits(height);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public int hashCode() {
        return cachedHashCode;
    }
}
