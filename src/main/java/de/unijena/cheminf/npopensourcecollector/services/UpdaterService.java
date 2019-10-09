package de.unijena.cheminf.npopensourcecollector.services;

import com.google.common.collect.Lists;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;

@Service
public class UpdaterService {

    @Autowired
    SourceNaturalProductRepository sourceNaturalProductRepository;

    @Autowired
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    List<Future<?>> futures = new ArrayList<Future<?>>();

    /*
    public void updateSourceNaturalProducts(){

        //get all sourceNaturalProduct
        System.out.println("Updating links...");

        List<SourceNaturalProduct> allSourceNaturalProducts = sourceNaturalProductRepository.findAll();

        for(SourceNaturalProduct snp : allSourceNaturalProducts){
            String unpid = snp.getUniqueNaturalProduct().getId();
            Optional<UniqueNaturalProduct> unp = uniqueNaturalProductRepository.findById(unpid);
            if(unp.isPresent()){
                UniqueNaturalProduct np = unp.get();
                snp.setUniqueNaturalProduct(np);
                sourceNaturalProductRepository.save(snp);
            }

        }
        System.out.println("done");
    }
*/


    public void updateSourceNaturalProductsParallelized(int nbThreads){

        System.out.println("Updating links...");

        //List<SourceNaturalProduct> allSourceNaturalProducts = sourceNaturalProductRepository.findAll();
        List<UniqueNaturalProduct> allUniqueNaturalProducts = uniqueNaturalProductRepository.findAll();


        try{

            System.out.println("Number of links to update: "+allUniqueNaturalProducts.size());

            ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(nbThreads);


            List<List<UniqueNaturalProduct>>  moleculeBatches =  Lists.partition(allUniqueNaturalProducts, 1000);

            List<Callable<Object>> todo = new ArrayList<Callable<Object>>(moleculeBatches.size());
            System.out.println("Total number of tasks:" + moleculeBatches.size());

            int taskcount = 0;



            for(List<UniqueNaturalProduct> oneBatch : moleculeBatches){

                UpdaterTask task = new UpdaterTask();

                task.setBatchOfMolecules(oneBatch);
                taskcount++;

                task.taskid=taskcount;

                Future<?> f = executor.submit(task);

                futures.add(f);

                //executor.execute(task);

                System.out.println("Task "+taskcount+" executing");

            }


            executor.shutdown();


        } catch (Exception e) {
        e.printStackTrace();
        }


        System.out.println("done");

    }


    public boolean processFinished(){

        boolean allFuturesDone = true;

        for(Future<?> future : this.futures){

            allFuturesDone &= future.isDone();

        }
        //System.out.println("Finished parallel computation of fragments");
        return allFuturesDone;
    }

}
