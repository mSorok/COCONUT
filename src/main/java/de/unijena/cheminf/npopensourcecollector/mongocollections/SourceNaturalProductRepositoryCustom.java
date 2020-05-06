package de.unijena.cheminf.npopensourcecollector.mongocollections;

import java.util.List;

public interface SourceNaturalProductRepositoryCustom {

    List<String> findUniqueInchiKeys();

    List<Object> findUniqueSourceNames();

}
