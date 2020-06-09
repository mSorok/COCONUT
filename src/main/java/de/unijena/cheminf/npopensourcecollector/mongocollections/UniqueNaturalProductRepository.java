package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.util.List;

public interface UniqueNaturalProductRepository  extends MongoRepository<UniqueNaturalProduct, String>, UniqueNaturalProductRepositoryCustom {

    public List<UniqueNaturalProduct> findBySmiles(String smiles);

    public List<UniqueNaturalProduct> findByInchi(String inchi);

    public List<UniqueNaturalProduct> findByInchikey(String inchikey);

    @Query("{ coconut_id : ?0}")
    public List<UniqueNaturalProduct> findByCoconut_id(String coconut_id);

    @Query("{ clean_smiles : ?0}")
    public List<UniqueNaturalProduct> findByClean_smiles(String clean_smiles);

    @Query("{molecular_formula : ?0}")
    public List<UniqueNaturalProduct> findByMolecular_formula(String molecular_formula);

    public List<UniqueNaturalProduct> findByName(String name);

    @Query("{ $text: { $search: ?0 } }")
    public List<UniqueNaturalProduct> fuzzyNameSearch(String name);




    @Query("{ npl_noh_score: { $exists:false } }")
    List<UniqueNaturalProduct> findAllByNPLScoreComputed();

    @Query("{ apol: { $exists:false } }")
    List<UniqueNaturalProduct> findAllByApolComputed();

    @Query("{ pubchemBits : { $bitsAllSet : ?0  }}")
    List<UniqueNaturalProduct> findAllPubchemBitsSet(byte[] querybits) ;




}
