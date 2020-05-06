package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.checkerframework.common.aliasing.qual.Unique;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.thymeleaf.util.StringUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Service
public class CreateCNPidService {

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;


    public static String prefix= "CNP";

    public void importIDs(String filename) {

        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader( new File(filename)));


            int count = 1;
            String line;


            while ((line = bufferedReader.readLine()) != null){
                ArrayList<String> dataline =new ArrayList<String>(Arrays.asList(line.split(","))); //coconut_id = 0, inchikey = 1

                List<UniqueNaturalProduct> unplist = uniqueNaturalProductRepository.findByInchikey(dataline.get(1));
                if(!unplist.isEmpty()){
                    for(UniqueNaturalProduct unp : unplist){
                        unp.setCoconut_id(dataline.get(0));
                    }

                }else{
                    System.out.println("BAD! Could not find "+dataline.get(1)+" in the new version of COCONUT!");
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

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

    public void createIDforNewMolecules(){

        int max_id = 0;
        List<UniqueNaturalProduct> allnp = uniqueNaturalProductRepository.findAll();

        ArrayList<UniqueNaturalProduct> unpWithoutId = new ArrayList<>();

        for(UniqueNaturalProduct np : allnp){

            if(np.getCoconut_id()==null || np.getCoconut_id()==""){
                unpWithoutId.add(np);

            }else{
                int coconut_tmp = Integer.parseInt( np.getCoconut_id().split("CNP")[1] );
                if(coconut_tmp>max_id ){
                    max_id = coconut_tmp;
                }
            }

        }


        max_id++;
        for(UniqueNaturalProduct ildnp : unpWithoutId){
            String coconut_id = prefix + StringUtils.repeat("0", 7-StringUtils.length(max_id)) + max_id;
            ildnp.setCoconut_id(coconut_id);
            uniqueNaturalProductRepository.save(ildnp);
            max_id++;
        }

    }
}
