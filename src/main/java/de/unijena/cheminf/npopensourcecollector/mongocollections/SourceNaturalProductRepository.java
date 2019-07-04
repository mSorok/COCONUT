package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.util.List;


public interface SourceNaturalProductRepository extends MongoRepository<SourceNaturalProduct, String>, SourceNaturalProductRepositoryCustom {

    List<SourceNaturalProduct> findBySimpleInchi(String inchi);

    List<SourceNaturalProduct> findBySimpleSmiles(String smiles);

    List<SourceNaturalProduct> findBySource(String source);

    List<SourceNaturalProduct> findBySimpleInchiKey(String inchikey);



}
