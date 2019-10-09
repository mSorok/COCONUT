package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.misc.BeanUtil;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.SourceNaturalProductRepository;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProductRepository;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.annotation.Transient;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Propagation;
import org.springframework.transaction.annotation.Transactional;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Optional;

@Service
@Transactional(propagation = Propagation.REQUIRED, readOnly = false)
public class UpdaterTask implements Runnable {

    @Autowired
    @Transient
    UniqueNaturalProductRepository uniqueNaturalProductRepository;

    @Autowired
    @Transient
    SourceNaturalProductRepository sourceNaturalProductRepository;

    List<UniqueNaturalProduct> batchOfMolecules;

    Integer taskid;


    @Override
    public void run() {
        this.uniqueNaturalProductRepository = BeanUtil.getBean(UniqueNaturalProductRepository.class);
        this.sourceNaturalProductRepository = BeanUtil.getBean(SourceNaturalProductRepository.class);


        for(UniqueNaturalProduct unp : batchOfMolecules){

            //find all SourceNaturalProducts that correspond to this UniqueNaturalproduct

            List<SourceNaturalProduct> allSNP = sourceNaturalProductRepository.findBySimpleInchiKey(unp.inchikey);
            for(SourceNaturalProduct snp : allSNP){
                snp.setUniqueNaturalProduct(unp);
                sourceNaturalProductRepository.save(snp);
            }

        }
        Date date = new Date();
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
        System.out.println("Task "+taskid+" finished at: "+formatter.format(date)+"\n");
    }


    public List<UniqueNaturalProduct> getBatchOfMolecules() {
        return batchOfMolecules;
    }

    public void setBatchOfMolecules(List<UniqueNaturalProduct> batchOfMolecules) {
        this.batchOfMolecules = batchOfMolecules;
    }

    public Integer getTaskid() {
        return taskid;
    }

    public void setTaskid(Integer taskid) {
        this.taskid = taskid;
    }
}
