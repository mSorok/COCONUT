package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.annotation.Id;
import org.springframework.data.mongodb.core.mapping.Document;

@Document
public class NPSimilarity {

    @Id
    public String id;

    public String uniqueNaturalProductID1;

    public String uniqueNaturalProductID2;

    public Double tanimoto;

    //public Double distanceMoment;


    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public Double getTanimoto() {
        return tanimoto;
    }

    public String getUniqueNaturalProductID1() {
        return uniqueNaturalProductID1;
    }

    public void setUniqueNaturalProductID1(String uniqueNaturalProductID1) {
        this.uniqueNaturalProductID1 = uniqueNaturalProductID1;
    }

    public String getUniqueNaturalProductID2() {
        return uniqueNaturalProductID2;
    }

    public void setUniqueNaturalProductID2(String uniqueNaturalProductID2) {
        this.uniqueNaturalProductID2 = uniqueNaturalProductID2;
    }

    public void setTanimoto(Double tanimoto) {
        this.tanimoto = tanimoto;
    }
/*
    public Double getDistanceMoment() {
        return distanceMoment;
    }

    public void setDistanceMoment(Double distanceMoment) {
        this.distanceMoment = distanceMoment;
    }
    */
}
