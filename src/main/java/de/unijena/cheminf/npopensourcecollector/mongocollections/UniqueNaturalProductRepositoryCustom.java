package de.unijena.cheminf.npopensourcecollector.mongocollections;

import java.util.List;

public interface UniqueNaturalProductRepositoryCustom {

    List<String> findAllCoconutIds();

    List<String> findAllInchiKeys();
}
