package de.unijena.cheminf.npopensourcecollector.mongocollections;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.domain.Sort;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.aggregation.Aggregation;
import org.springframework.data.mongodb.core.aggregation.AggregationResults;
import org.springframework.data.mongodb.core.aggregation.GroupOperation;
import org.springframework.data.mongodb.core.query.Criteria;

import java.util.List;

import static org.springframework.data.mongodb.core.aggregation.Aggregation.*;

public class SourceNaturalProductRepositoryImpl implements SourceNaturalProductRepositoryCustom {

    private final MongoTemplate mongoTemplate;

    @Autowired
    public SourceNaturalProductRepositoryImpl(MongoTemplate mongoTemplate) {
        this.mongoTemplate = mongoTemplate;
    }


    @Override
    public List<String> findUniqueInchiKeys(){

        GroupOperation groupByInchikey = group("simpleInchiKey");

        Aggregation aggregation = newAggregation(groupByInchikey);
        AggregationResults<String> groupResults = mongoTemplate.aggregate(aggregation, "sourceNaturalProduct", String.class);

        List<String> result =  groupResults.getMappedResults();

        return result;
    }
    /*public List<Object> findUniqueInchiKeys() {
        return mongoTemplate.query(SourceNaturalProduct.class).distinct("simpleInchiKey").all()  ;
    }*/

    @Override
    public List<Object> findUniqueSourceNames(){
        return mongoTemplate.query(SourceNaturalProduct.class).distinct("source").all() ;
    }


}
