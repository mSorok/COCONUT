package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.mongodb.repository.MongoRepository;

import java.util.List;

public interface UniqueNaturalProductRepository  extends MongoRepository<UniqueNaturalProduct, String> {

    public List<UniqueNaturalProduct> findBySmiles(String smiles);

    public List<UniqueNaturalProduct> findByInchi(String inchi);

    public List<UniqueNaturalProduct> findByInchikey(String inchikey);
}
