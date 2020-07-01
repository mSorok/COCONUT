package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.annotation.Id;
import org.springframework.data.mongodb.repository.MongoRepository;

public interface PubFingerprintsCountsRepository extends MongoRepository<PubFingerprintsCounts, String> {



}
