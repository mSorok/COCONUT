package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.util.List;


public interface SourceNaturalProductRepository extends MongoRepository<SourceNaturalProduct, String>, SourceNaturalProductRepositoryCustom {


    List<SourceNaturalProduct> findBySimpleInchiKey(String inchikey);




}
