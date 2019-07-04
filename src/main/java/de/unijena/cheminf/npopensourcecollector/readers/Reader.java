package de.unijena.cheminf.npopensourcecollector.readers;


import org.javatuples.Pair;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.File;
import java.util.ArrayList;

public interface Reader {

    void readFile(File file);

    ArrayList<IAtomContainer> returnCorrectMolecules();

    String returnSource();
}
