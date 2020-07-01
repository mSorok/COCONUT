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


    public void clearIDs(){

        List<UniqueNaturalProduct> allunp = uniqueNaturalProductRepository.findAll();

        for(UniqueNaturalProduct unp : allunp){

            unp.setCoconut_id("");

            uniqueNaturalProductRepository.save(unp);


        }

    }

    public void importIDs(String filename) {

        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader( new File(filename)));


            int count = 1;
            String line;


            while ((line = bufferedReader.readLine()) != null){
                if(!line.startsWith("coconut_id")) {
                    //ArrayList<String> dataline = new ArrayList<String>(Arrays.asList(line.split(","))); //coconut_id = 0, inchikey = 1
                    String[] dataTab = line.split(",") ;

                    List<UniqueNaturalProduct> unplist = uniqueNaturalProductRepository.findByInchikey(dataTab[1]);


                    if (!unplist.isEmpty()) {
                        for (UniqueNaturalProduct unp : unplist) {
                            unp.setCoconut_id(dataTab[0]);
                            uniqueNaturalProductRepository.save(unp);
                        }

                    } else {
                        System.out.println("BAD! Could not find " + dataTab[0] + " in the new version of COCONUT!");
                    }


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

            if(np.coconut_id.equals("")){
                unpWithoutId.add(np);

            }else if(np.coconut_id.startsWith("CNP")){
                int coconut_tmp = Integer.parseInt( np.getCoconut_id().split("NP")[1] );
                if(coconut_tmp>max_id ){
                    max_id = coconut_tmp;
                }
            }

        }


        max_id+=1;
        for(UniqueNaturalProduct ildnp : unpWithoutId){
            String coconut_id = prefix + StringUtils.repeat("0", 7-StringUtils.length(max_id)) + max_id;
            ildnp.setCoconut_id(coconut_id);
            uniqueNaturalProductRepository.save(ildnp);
            max_id+=1;
        }

    }
}
