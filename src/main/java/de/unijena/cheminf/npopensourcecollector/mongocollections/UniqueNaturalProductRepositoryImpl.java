package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.mongodb.core.MongoTemplate;

import java.util.List;

public class UniqueNaturalProductRepositoryImpl implements UniqueNaturalProductRepositoryCustom {

    private final MongoTemplate mongoTemplate;

    @Autowired
    public UniqueNaturalProductRepositoryImpl(MongoTemplate mongoTemplate) {
        this.mongoTemplate = mongoTemplate;
    }


    @Override
    public List<String> findAllCoconutIds() {

        List<String> coconut_ids_list = mongoTemplate.query(UniqueNaturalProduct.class)
                .distinct("coconut_id")
                .as(String.class)
                .all();
        return coconut_ids_list;
    }
}
