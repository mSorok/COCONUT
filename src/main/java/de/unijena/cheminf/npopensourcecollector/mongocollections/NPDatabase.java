package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.annotation.Id;
import org.springframework.data.mongodb.core.mapping.Document;

@Document
public class NPDatabase {

    @Id
    public String id;

    String name;

    String localFileName;

    String url;

    String comments;

    Integer nb_unique_molecules;


    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
    }

    public String getComments() {
        return comments;
    }

    public void setComments(String comments) {
        this.comments = comments;
    }

    public Integer getNb_unique_molecules() {
        return nb_unique_molecules;
    }

    public void setNb_unique_molecules(Integer nb_unique_molecules) {
        this.nb_unique_molecules = nb_unique_molecules;
    }

    public String getLocalFileName() {
        return localFileName;
    }

    public void setLocalFileName(String localFileName) {
        this.localFileName = localFileName;
    }
}
