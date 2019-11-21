package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.checkerframework.common.aliasing.qual.Unique;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.thymeleaf.util.StringUtils;

import java.util.ArrayList;
import java.util.List;

@Service
public class CreateCNPidService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;


    public static String prefix= "CNP";

    public void createDeNovoIDs(){

        List<UniqueNaturalProduct> allunp = uniqueNaturalProductRepository.findAll();

        int count = 1;
        for(UniqueNaturalProduct unp : allunp){


            String coconut_id = prefix + StringUtils.repeat("0", 7-StringUtils.length(count)) + count;

            unp.setCoconut_id(coconut_id);

            uniqueNaturalProductRepository.save(unp);

            count++;

        }

    }

    public void addID(){

        //TODO

    }
}
