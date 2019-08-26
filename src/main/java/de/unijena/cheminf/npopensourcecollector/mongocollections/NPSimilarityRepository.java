package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.data.mongodb.repository.MongoRepository;

public interface NPSimilarityRepository extends MongoRepository<NPSimilarity, String> {
}
