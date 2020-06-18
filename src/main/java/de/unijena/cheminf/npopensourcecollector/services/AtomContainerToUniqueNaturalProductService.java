package de.unijena.cheminf.npopensourcecollector.services;

import de.unijena.cheminf.npopensourcecollector.mongocollections.UniqueNaturalProduct;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.springframework.stereotype.Service;

@Service
public class AtomContainerToUniqueNaturalProductService {

    public IAtomContainer createAtomContainer(UniqueNaturalProduct unp){

        IAtomContainer ac = null;

        try {
            SmilesParser sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
            ac   = sp.parseSmiles( unp.smiles );

            if(ac == null || ac.isEmpty()){
                System.out.println("Failed to recreate AC for NP "+unp.inchikey);
            }
        } catch (InvalidSmilesException e) {
            System.err.println(e.getMessage());
        }

        return ac;

    }
}
