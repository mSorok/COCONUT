package de.unijena.cheminf.npopensourcecollector.mongocollections;

import java.util.ArrayList;

public class UncomplicatedTaxonomy {

    String cleaned_organism_id;

    String organism_value;

    String organism_url;

    public String wikidata_id;


    public String reference_wikidata_id;


    String domain = null;

    String superkingdom = null;

    String kingdom = null;

    String phylum = null;

    String classx = null;

    String order = null;

    String family = null;

    String genus = null;

    String species = null;


    public UncomplicatedTaxonomy() {
    }




    public String getAllRanks(){
        ArrayList taxRanks = new ArrayList();
        String prettyString = "";



        if(this.domain != null){
            taxRanks.add("domain");
        }

        if(this.superkingdom != null){
            taxRanks.add("superkingdom");
        }

        if(this.kingdom != null){
            taxRanks.add("kingdom");
        }

        if (this.phylum != null){
            taxRanks.add("phylum");
        }

        if(this.classx != null){
            taxRanks.add("classx");
        }

        if(this.order != null){
            taxRanks.add("order");
        }

        if(this.family != null){
            taxRanks.add("family");
        }

        if(this.genus != null){
            taxRanks.add("genus");
        }

        if(this.species != null){
            taxRanks.add("species");
        }


        prettyString = String.join(" | ", taxRanks);
        return prettyString;
    }



    @Override
    public String toString(){

        ArrayList taxNames = new ArrayList();
        String prettyString = "";

        if(this.domain != null){
            taxNames.add(this.domain);
        }

        if(this.superkingdom != null){
            taxNames.add(this.superkingdom);
        }

        if(this.kingdom != null){
            taxNames.add(this.kingdom);
        }

        if (this.phylum != null){
            taxNames.add(this.phylum);
        }

        if(this.classx != null){
            taxNames.add(this.classx);
        }

        if(this.order != null){
            taxNames.add(this.order);
        }

        if(this.family != null){
            taxNames.add(this.family);
        }

        if(this.genus != null){
            taxNames.add(this.genus);
        }

        if(this.species != null){
            taxNames.add(this.species);
        }



        prettyString = String.join(" | ", taxNames);
        return prettyString;
    }



    public String getWikidata_id() {
        return wikidata_id;
    }

    public void setWikidata_id(String wikidata_id) {
        this.wikidata_id = wikidata_id;
    }

    public String getReference_wikidata_id() {
        return reference_wikidata_id;
    }

    public void setReference_wikidata_id(String reference_wikidata_id) {
        this.reference_wikidata_id = reference_wikidata_id;
    }

    public String getCleaned_organism_id() {
        return cleaned_organism_id;
    }

    public void setCleaned_organism_id(String cleaned_organism_id) {
        this.cleaned_organism_id = cleaned_organism_id;
    }

    public String getOrganism_value() {
        return organism_value;
    }

    public void setOrganism_value(String organism_value) {
        this.organism_value = organism_value;
    }

    public String getOrganism_url() {
        return organism_url;
    }

    public void setOrganism_url(String organism_url) {
        this.organism_url = organism_url;
    }

    public String getDomain() {
        return domain;
    }

    public void setDomain(String domain) {
        this.domain = domain;
    }

    public String getSuperkingdom() {
        return superkingdom;
    }

    public void setSuperkingdom(String superkingdom) {
        this.superkingdom = superkingdom;
    }

    public String getKingdom() {
        return kingdom;
    }

    public void setKingdom(String kingdom) {
        this.kingdom = kingdom;
    }

    public String getPhylum() {
        return phylum;
    }

    public void setPhylum(String phylum) {
        this.phylum = phylum;
    }

    public String getClassx() {
        return classx;
    }

    public void setClassx(String classx) {
        this.classx = classx;
    }

    public String getOrder() {
        return order;
    }

    public void setOrder(String order) {
        this.order = order;
    }

    public String getFamily() {
        return family;
    }

    public void setFamily(String family) {
        this.family = family;
    }

    public String getGenus() {
        return genus;
    }

    public void setGenus(String genus) {
        this.genus = genus;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }
}
