package de.unijena.cheminf.npopensourcecollector.misc;

import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Class to discard counter ions and other small disconnected fragments
 *
 * @author kalai
 */

@Service
public class MoleculeConnectivityChecker {

    public List<IAtomContainer> checkConnectivity(IAtomContainer molecule) {

        List<IAtomContainer> curated = new ArrayList<IAtomContainer>();
        Map<Object, Object> properties = molecule.getProperties();
        IAtomContainerSet mols = ConnectivityChecker.partitionIntoMolecules(molecule);
        for (int i = 0; i < mols.getAtomContainerCount(); i++) {
            mols.getAtomContainer(i).setProperties(properties);
            curated.add(mols.getAtomContainer(i));
        }
        return curated;
    }
}
