package objects;

/**
 *
 * @author devita
 */
public enum CVTType {
    UDF(0, "UDF"),
    SVT(1, "SVT"),
    BMTC(2, "BMTC"),
    BMTZ(3, "BMTZ");

    private final int id;
    private final String name;

    CVTType() {
        id = 0;
        name = "UDF";
    }

    CVTType(int id, String name) {
        this.id = id;
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public int getId() {
        return id;
    }
}