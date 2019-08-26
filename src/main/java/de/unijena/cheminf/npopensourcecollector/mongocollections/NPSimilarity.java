package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.annotation.Id;
import org.springframework.data.mongodb.core.mapping.Document;

@Document
public class NPSimilarity {

    @Id
    public String id;

    public UniqueNaturalProduct uniqueNaturalProduct1;

    public UniqueNaturalProduct uniqueNaturalProduct2;

    public Double tanimoto;

    //public Double distanceMoment;


    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public UniqueNaturalProduct getUniqueNaturalProduct1() {
        return uniqueNaturalProduct1;
    }

    public void setUniqueNaturalProduct1(UniqueNaturalProduct uniqueNaturalProduct1) {
        this.uniqueNaturalProduct1 = uniqueNaturalProduct1;
    }

    public UniqueNaturalProduct getUniqueNaturalProduct2() {
        return uniqueNaturalProduct2;
    }

    public void setUniqueNaturalProduct2(UniqueNaturalProduct uniqueNaturalProduct2) {
        this.uniqueNaturalProduct2 = uniqueNaturalProduct2;
    }

    public Double getTanimoto() {
        return tanimoto;
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
