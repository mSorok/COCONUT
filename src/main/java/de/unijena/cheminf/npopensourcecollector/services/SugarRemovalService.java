package de.unijena.cheminf.npopensourcecollector.services;




/**
 * @author Jonas Schaub
 * @author Maria Sorokina
 */


import de.unijena.cheminf.npopensourcecollector.misc.MoleculeConnectivityChecker;
import de.unijena.cheminf.npopensourcecollector.readers.ReaderService;
import net.sf.jniinchi.INCHI_OPTION;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerComparator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;


import javax.swing.plaf.synth.SynthEditorPaneUI;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

@Service
public class SugarRemovalService {


    public static enum StructuresToKeepMode {
        ALL (0),
        HEAVY_ATOM_COUNT (5),
        MOLECULAR_WEIGHT (60);
        private final int defaultThreshold;
        StructuresToKeepMode(int aDefaultValue) {
            this.defaultThreshold = aDefaultValue;
        }
        public int getDefaultThreshold() {
            return this.defaultThreshold;
        }
    }



    public static final String[] LINEAR_SUGARS_SMILES = {
            //*aldoses*
            //note: no octose and so on
            "C(C(C(C(C(C(C=O)O)O)O)O)O)O", //aldoheptose
            "C(C(C(C(C(C=O)O)O)O)O)O", //aldohexose
            "C(C(C(C(C=O)O)O)O)O", //aldopentose
            "C(C(C(C=O)O)O)O", //aldotetrose
            "C(C(C=O)O)O", //aldotriose
            //*ketoses*
            //note: no octose and so on
            "C(C(C(C(C(C(CO)O)O)O)O)=O)O", //2-ketoheptose
            "C(C(C(C(C(CO)O)O)O)=O)O", //2-ketohexose
            "C(C(C(C(CO)O)O)=O)O", //2-ketopentose
            "C(C(C(CO)O)=O)O", //2-ketotetrose
            "C(C(CO)=O)O", //2-ketotriose
            //*sugar alcohols*
            //note: no octitol and so on
            "C(C(C(C(C(C(CO)O)O)O)O)O)O", //heptitol
            "C(C(C(C(C(CO)O)O)O)O)O", //hexitol
            "C(C(C(C(CO)O)O)O)O", //pentitol
            "C(C(C(CO)O)O)O", //tetraitol
            "C(C(CO)O)O", //triol
            //*deoxy sugars*
            "C(C(C(C(CC=O)O)O)O)O" //2-deoxyhexose
    };

    public static final String [] RING_SUGARS_SMILES = {
            "C1CCOC1", //tetrahydrofuran to match all 5-membered sugar rings
            "C1CCOCC1", //tetrahydropyran to match all 6-membered sugar rings
            "C1CCCOCC1" //oxepane to match all 7-membered sugar rings
    };

    public static final boolean REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT = false;
    public static final boolean DETECT_GLYCOSIDIC_BOND_DEFAULT = false;
    public static final boolean REMOVE_ONLY_TERMINAL_DEFAULT = true;
    public static final boolean INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT = true;
    public static final double ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT = 0.5;
    public static final boolean SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT = true;
    public static final String INDEX_PROPERTY_KEY = "SUGAR_REMOVAL_UTILITY_INDEX";
    public static final String CONTAINS_LINEAR_SUGAR_PROPERTY_KEY = "CONTAINS_LINEAR_SUGAR";
    public static final String CONTAINS_SUGAR_PROPERTY_KEY = "CONTAINS_SUGAR";
    public static final String CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY = "CONTAINS_CIRCULAR_SUGAR";
    public static final int LINEAR_SUGAR_CANDIDATE_MIN_SIZE_DEFAULT = 4;
    public static final int LINEAR_SUGAR_CANDIDATE_MAX_SIZE_DEFAULT = 7;
    public static final SmartsPattern ESTER_SMARTS_PATTERN = SmartsPattern.create("[C](=O)-[O!R]-[C]");
    public static final SmartsPattern ETHER_SMARTS_PATTERN = SmartsPattern.create("[C]-[O!R]-[C]");
    public static final SmartsPattern PEROXIDE_SMARTS_PATTERN = SmartsPattern.create("[C]-[O!R]-[O!R]-[C]");


    UniversalIsomorphismTester universalIsomorphismTester ;

    @Autowired
    MoleculeConnectivityChecker mcc;

    @Autowired
    ReaderService readerService;




    private List<IAtomContainer> ringSugars;
    private List<IAtomContainer> linearSugars;

    private List<DfPattern> linearSugarPatterns;

    private StructuresToKeepMode structuresToKeepMode;
    public static final StructuresToKeepMode STRUCTURES_TO_KEEP_MODE_DEFAULT = StructuresToKeepMode.HEAVY_ATOM_COUNT;
    private int structureToKeepModeThreshold;


    private boolean detectGlycosidicBond;
    private boolean removeOnlyTerminal;
    private boolean removeLinearSugarsInRing;
    private boolean includeNrOfAttachedOxygens;
    private double attachedOxygensToAtomsInRingRatioThreshold;
    private int linearSugarCandidateMinSize;
    private int linearSugarCandidateMaxSize;


    private boolean setPropertyOfSugarContainingMolecules;




    IAtomContainer removeSugarsFromAtomContainer(IAtomContainer  moleculeToProcess){


        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);

        List<IAtomContainer> results = null ;



        List options = new ArrayList();
        options.add(INCHI_OPTION.SNon);
        options.add(INCHI_OPTION.ChiralFlagOFF);
        options.add(INCHI_OPTION.AuxNone);




        //removing all sugars
        try {

            setRemoveLinearSugarsInRing(false);
            setPropertyOfSugarContainingMolecules(true);
            setRemoveOnlyTerminalSugars(false);


            try {
                results = removeAndReturnCircularAndLinearSugars(moleculeToProcess, false);

                moleculeToProcess = results.get(0);
                //the molecule to process can be in several parts: need to separate them
                //List<IAtomContainer> listAC = mcc.checkConnectivity(moleculeToProcess);

                //returning only the biggest fragment
                moleculeToProcess = selectBiggestUnconnectedFragment(moleculeToProcess);


            } catch (CloneNotSupportedException  e) {
                e.printStackTrace();
                return null;
            }

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }

        //add a nice way of dealing with long smiles



        return moleculeToProcess;
    }





    /******************** Actual sugar removal methods *****************************/

    protected void prepareSugars(){
        universalIsomorphismTester = new UniversalIsomorphismTester();

        this.linearSugars = new ArrayList<>(LINEAR_SUGARS_SMILES.length);
        this.ringSugars = new ArrayList<>(RING_SUGARS_SMILES.length);
        this.linearSugarPatterns = new ArrayList<>(LINEAR_SUGARS_SMILES.length);

        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        //adding linear sugars to list
        for (String tmpSmiles : LINEAR_SUGARS_SMILES) {
            try {
                this.linearSugars.add(tmpSmilesParser.parseSmiles(tmpSmiles));
            } catch (Exception anException) {
                anException.printStackTrace();
            }
        }
        //sorting for size decreasing; the patterns parsed afterwards are not sorted that easily, so sorting is done now
        Comparator<IAtomContainer> tmpComparator = new AtomContainerComparator().reversed();
        //note: this can throw various exceptions but they should not appear here
        this.linearSugars.sort(tmpComparator);
        //adding ring sugars to list
        for (String tmpSmiles : RING_SUGARS_SMILES) {
            try {
                this.ringSugars.add(tmpSmilesParser.parseSmiles(tmpSmiles));
            } catch (Exception anException) {
                anException.printStackTrace();
            }
        }
        this.ringSugars.sort(tmpComparator);

        //parsing linear sugars into patterns
        for(IAtomContainer tmpSugarAC : this.linearSugars){
            try {
                this.linearSugarPatterns.add(DfPattern.findSubstructure(tmpSugarAC));
            } catch (Exception anException) {
                anException.printStackTrace();
            }
        }



        this.detectGlycosidicBond = DETECT_GLYCOSIDIC_BOND_DEFAULT;
        this.removeOnlyTerminal = REMOVE_ONLY_TERMINAL_DEFAULT;
        this.structuresToKeepMode = STRUCTURES_TO_KEEP_MODE_DEFAULT;
        this.structureToKeepModeThreshold = this.structuresToKeepMode.defaultThreshold;
        this.includeNrOfAttachedOxygens = INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT;
        this.attachedOxygensToAtomsInRingRatioThreshold =
                ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT;
        this.removeLinearSugarsInRing = REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT;
        this.setPropertyOfSugarContainingMolecules = SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT;
        this.linearSugarCandidateMinSize = LINEAR_SUGAR_CANDIDATE_MIN_SIZE_DEFAULT;
        this.linearSugarCandidateMaxSize = LINEAR_SUGAR_CANDIDATE_MAX_SIZE_DEFAULT;

    }


    public List<String> getLinearSugars() {
        List<String> tmpSmilesList = new ArrayList<>(this.linearSugars.size());
        SmilesGenerator tmpSmilesGen = new SmilesGenerator(SmiFlavor.Unique);
        for (IAtomContainer tmpLinearSugar : this.linearSugars) {
            String tmpSmiles = null;
            try {
                tmpSmiles = tmpSmilesGen.create(tmpLinearSugar);
            } catch (CDKException aCDKException) {
                aCDKException.printStackTrace();
            }
            if (!Objects.isNull(tmpSmiles)) {
                try {
                    tmpSmilesList.add(tmpSmiles);
                } catch (Exception anException) {
                    anException.printStackTrace();
                }
            }
        }
        return tmpSmilesList;
    }

    public List<String> getCircularSugars() {
        List<String> tmpSmilesList = new ArrayList<>(this.ringSugars.size());
        SmilesGenerator tmpSmilesGen = new SmilesGenerator(SmiFlavor.Unique);
        for (IAtomContainer tmpRingSugar : this.ringSugars) {
            String tmpSmiles = null;
            try {
                tmpSmiles = tmpSmilesGen.create(tmpRingSugar);
            } catch (CDKException aCDKException) {
                aCDKException.printStackTrace();
            }
            if (!Objects.isNull(tmpSmiles)) {
                try {
                    tmpSmilesList.add(tmpSmiles);
                } catch (Exception anException) {
                    anException.printStackTrace();
                }
            }
        }
        return tmpSmilesList;
    }

    public boolean isGlycosidicBondDetected() {
        return this.detectGlycosidicBond;
    }

    public boolean areOnlyTerminalSugarsRemoved() {
        return this.removeOnlyTerminal;
    }

    public StructuresToKeepMode getStructuresToKeepMode() {
        return this.structuresToKeepMode;
    }

    public int getStructureToKeepModeThreshold() {
        return this.structureToKeepModeThreshold;
    }

    public boolean isNrOfAttachedOxygensIncluded() {
        return this.includeNrOfAttachedOxygens;
    }

    public double getAttachedOxygensToAtomsInRingRatioThreshold() {
        return this.attachedOxygensToAtomsInRingRatioThreshold;
    }

    public boolean areLinearSugarsInRingsRemoved() {
        return this.removeLinearSugarsInRing;
    }

    public boolean isPropertyOfSugarContainingMoleculesSet() {
        return this.setPropertyOfSugarContainingMolecules;
    }

    public boolean arePropertiesOfSugarContainingMoleculesSet() {
        return this.setPropertyOfSugarContainingMolecules;
    }



    public int getLinearSugarCandidateMinSize() {
        return this.linearSugarCandidateMinSize;
    }


    public int getLinearSugarCandidateMaxSize() {
        return this.linearSugarCandidateMaxSize;
    }

    public void clearCircularSugars() {
        try {
            this.ringSugars.clear();
        } catch (UnsupportedOperationException anException) {
            anException.printStackTrace();
            this.ringSugars = new ArrayList<>(RING_SUGARS_SMILES.length);
        }
    }

    public void clearLinearSugars() {
        try {
            this.linearSugars.clear();
            this.linearSugarPatterns.clear();
        } catch (UnsupportedOperationException anException) {
            anException.printStackTrace();
            this.linearSugars = new ArrayList<>(LINEAR_SUGARS_SMILES.length);
            this.linearSugarPatterns = new ArrayList<>(LINEAR_SUGARS_SMILES.length);
        }
    }

    public void setDetectGlycosidicBond(boolean aBoolean) {
        this.detectGlycosidicBond = aBoolean;
    }

    public void setRemoveOnlyTerminalSugars(boolean aBoolean) {
        this.removeOnlyTerminal = aBoolean;
    }

    public void setIncludeNrOfAttachedOxygens(boolean aBoolean) {
        this.includeNrOfAttachedOxygens = aBoolean;
    }

    public void setRemoveLinearSugarsInRing(boolean aBoolean) {
        this.removeLinearSugarsInRing = aBoolean;
    }

    public void setPropertyOfSugarContainingMolecules(boolean aBoolean) {
        this.setPropertyOfSugarContainingMolecules = aBoolean;
    }

    public void setLinearSugarCandidateMinSize(int aMinSize) throws IllegalArgumentException {
        if (aMinSize < 1) {
            throw new IllegalArgumentException("Given minimum size is smaller than 1.");
        }
        this.linearSugarCandidateMinSize = aMinSize;
    }

    public void setLinearSugarCandidateMaxSize(int aMaxSize) throws IllegalArgumentException {
        if (aMaxSize < 1) {
            throw new IllegalArgumentException("Given maximum size is smaller than 1.");
        }
        this.linearSugarCandidateMaxSize = aMaxSize;
    }

    public void restoreDefaultSettings() {
        this.detectGlycosidicBond = DETECT_GLYCOSIDIC_BOND_DEFAULT;
        this.removeOnlyTerminal = REMOVE_ONLY_TERMINAL_DEFAULT;
        this.structuresToKeepMode = STRUCTURES_TO_KEEP_MODE_DEFAULT;
        this.structureToKeepModeThreshold = this.structuresToKeepMode.defaultThreshold;
        this.includeNrOfAttachedOxygens = INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT;
        this.attachedOxygensToAtomsInRingRatioThreshold =
                ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT;
        this.removeLinearSugarsInRing = REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT;
        this.setPropertyOfSugarContainingMolecules = SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT;
        this.linearSugarCandidateMinSize = LINEAR_SUGAR_CANDIDATE_MIN_SIZE_DEFAULT;
        this.linearSugarCandidateMaxSize = LINEAR_SUGAR_CANDIDATE_MAX_SIZE_DEFAULT;
    }


    public boolean hasLinearSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        //throws NullPointerException if molecule is null
        this.setIndices(aMolecule);
        //throws NullPointerException if molecule is null
        List<IAtomContainer> tmpSugarCandidates = this.getLinearSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            aMolecule.setProperty(CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        return tmpContainsSugar;
    }

    public boolean hasCircularSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        //throws NullPointerException if molecule is null
        this.setIndices(aMolecule);
        //throws NullPointerException if molecule is null
        List<IAtomContainer> tmpSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            aMolecule.setProperty(CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        return tmpContainsSugar;
    }


    public boolean hasSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        //throws NullPointerException if molecule is null
        this.setIndices(aMolecule);
        //throws NullPointerException if molecule is null
        List<IAtomContainer> tmpCircularSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsCircularSugar = !tmpCircularSugarCandidates.isEmpty();
        //throws NullPointerException if molecule is null
        List<IAtomContainer> tmpLinearSugarCandidates = this.getLinearSugarCandidates(aMolecule);
        boolean tmpContainsLinearSugar = !tmpLinearSugarCandidates.isEmpty();
        boolean tmpContainsSugar = (tmpContainsCircularSugar || tmpContainsLinearSugar);
        if (this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            aMolecule.setProperty(CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsCircularSugar);
            aMolecule.setProperty(CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsLinearSugar);
        }
        return tmpContainsSugar;
    }

    public int getNumberOfCircularSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return 0;
        }
        List<IAtomContainer> tmpCircularSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        int tmpSize = tmpCircularSugarCandidates.size();
        return tmpSize;
    }


    public IAtomContainer removeCircularSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        IAtomContainer tmpNewMolecule = this.removeAndReturnCircularSugars(aMolecule, aShouldBeCloned).get(0);
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }


    public boolean removeCircularSugars(IAtomContainer aMolecule)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        List<IAtomContainer> tmpDeglycosylatedMoleculeAndSugarMoietiesList =
                this.removeAndReturnCircularSugars(aMolecule, false);
        boolean tmpSomethingWasRemoved = (tmpDeglycosylatedMoleculeAndSugarMoietiesList.size() > 1);
        return tmpSomethingWasRemoved;
    }


    public List<IAtomContainer> removeAndReturnCircularSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        //throws NullPointerException if molecule is null
        this.setIndices(tmpNewMolecule);
        //throws NullPointerException if molecule is null
        List<IAtomContainer> tmpSugarCandidates = this.getCircularSugarCandidates(tmpNewMolecule);
        /*note: this means that there are matches of the circular sugar patterns and that they adhere to most of
        the given settings. The exception is that they might not be terminal*/
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            tmpNewMolecule.setProperty(CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            tmpNewMolecule.setProperty(CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        List<IAtomContainer> tmpResultList = new ArrayList<>(tmpSugarCandidates.size() + 1);
        tmpResultList.add(0, tmpNewMolecule);
        if (tmpContainsSugar) {
            //throws NullPointerException and IllegalArgumentException
            tmpResultList.addAll(1, this.removeSugarCandidates(tmpNewMolecule, tmpSugarCandidates));
        }
        //the molecule at index 0 may be empty and may be unconnected, based on the settings
        return tmpResultList;
    }




    public List<IAtomContainer> removeSugarCandidates(IAtomContainer aMolecule, List<IAtomContainer> aCandidateList)
            throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty() || aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        for (IAtomContainer tmpSubstructure : aCandidateList) {
            boolean tmpIsParent = true;
            if (Objects.isNull(tmpSubstructure) || tmpSubstructure.isEmpty()) {
                continue;
            }
            for (IAtom tmpAtom : tmpSubstructure.atoms()) {
                if (!aMolecule.contains(tmpAtom)) {
                    tmpIsParent = false;
                    break;
                }
            }
            if (!tmpIsParent) {
                throw new IllegalArgumentException("At least one of the possible sugar-like substructures is not " +
                        "actually part of the given molecule.");
            }
        }
        // a copy of the list is needed to avoid iterating over the same elements again if only terminal moieties are removed
        List<IAtomContainer> tmpSugarCandidates = new ArrayList(aCandidateList);
        // the to be returned list of removed moieties
        List<IAtomContainer> tmpRemovedSugarMoieties = new ArrayList<>(aCandidateList.size());
        if (this.removeOnlyTerminal) {
            //Only terminal sugars should be removed
            //but the definition of terminal depends on the set structures to keep mode!
            //decisions based on this setting are made in the respective private method
            //No unconnected structures result at the end or at an intermediate step
            boolean tmpContainsNoTerminalSugar = false;
            while (!tmpContainsNoTerminalSugar) {
                boolean tmpSomethingWasRemoved = false;
                for (int i = 0; i < tmpSugarCandidates.size(); i++) {
                    IAtomContainer tmpCandidate = tmpSugarCandidates.get(i);
                    if (Objects.isNull(tmpCandidate) || tmpCandidate.isEmpty()) {
                        continue;
                    }
                    boolean tmpIsTerminal = false;
                    try {
                        //also throws NullPointerExceptions or IllegalArgumentExceptions but they are simply passed on
                        // by this calling method
                        tmpIsTerminal = this.isTerminal(tmpCandidate, aMolecule, tmpSugarCandidates);
                    } catch (CloneNotSupportedException aCloneNotSupportedException) {
                        aCloneNotSupportedException.printStackTrace();

                        throw new IllegalArgumentException("Could not clone one candidate sugar structure and therefore " +
                                "not determine whether it is terminal or not.");
                    }
                    if (tmpIsTerminal) {
                        for (IAtom tmpAtom : tmpCandidate.atoms()) {
                            if (aMolecule.contains(tmpAtom)) {
                                aMolecule.removeAtom(tmpAtom);
                            }
                        }
                        tmpRemovedSugarMoieties.add(tmpCandidate);
                        tmpSugarCandidates.remove(i);
                        //The removal shifts the remaining indices!
                        i = i - 1;
                        if (!aMolecule.isEmpty()) {
                            //to clear away leftover unconnected fragments that are not to be kept due to the settings and
                            // to generate valid valences by adding implicit hydrogen atoms
                            //throws NullPointerException if molecule is null
                            this.postProcessAfterRemoval(aMolecule);
                        }
                        //atom container may be empty after that
                        if (aMolecule.isEmpty()) {
                            tmpContainsNoTerminalSugar = true;
                            break;
                        }
                        tmpSomethingWasRemoved = true;
                    }
                }
                if (!tmpSomethingWasRemoved) {
                    tmpContainsNoTerminalSugar = true;
                }
            }
        } else {
            //all sugar moieties are removed, may result in an unconnected atom container
            for (IAtomContainer tmpSugarCandidate : tmpSugarCandidates) {
                for (IAtom tmpAtom : tmpSugarCandidate.atoms()) {
                    if (aMolecule.contains(tmpAtom)) {
                        aMolecule.removeAtom(tmpAtom);
                    }
                }
                tmpRemovedSugarMoieties.add(tmpSugarCandidate);
            }
        }
        if (!aMolecule.isEmpty()) {
            //to clear away leftover unconnected fragments that are not to be kept due to the settings and
            // to generate valid valences by adding implicit hydrogen atoms
            //throws NullPointerException if molecule is null
            this.postProcessAfterRemoval(aMolecule);
        }
        return tmpRemovedSugarMoieties;
    }





    public IAtomContainer removeLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        IAtomContainer tmpNewMolecule = this.removeAndReturnLinearSugars(aMolecule, aShouldBeCloned).get(0);
        //the molecule at index 0 may be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }


    public boolean removeLinearSugars(IAtomContainer aMolecule)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        List<IAtomContainer> tmpDeglycosylatedMoleculeAndSugarMoietiesList =
                this.removeAndReturnLinearSugars(aMolecule, false);
        boolean tmpSomethingWasRemoved = (tmpDeglycosylatedMoleculeAndSugarMoietiesList.size() > 1);
        return tmpSomethingWasRemoved;
    }


    public List<IAtomContainer> removeAndReturnLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        //</editor-fold>
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        //throws NullPointerException if molecule is null
        this.setIndices(tmpNewMolecule);
        //throws NullPointerException if molecule is null
        List<IAtomContainer> tmpSugarCandidates = this.getLinearSugarCandidates(tmpNewMolecule);
        /*note: this means that there are matches of the linear sugar patterns and that they adhere to most of
        the given settings. The exception is that they might not be terminal*/
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            tmpNewMolecule.setProperty(CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            tmpNewMolecule.setProperty(CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        List<IAtomContainer> tmpResultList = new ArrayList<>(tmpSugarCandidates.size() + 1);
        tmpResultList.add(0, tmpNewMolecule);
        if (tmpContainsSugar) {
            //throws NullPointerException and IllegalArgumentException
            tmpResultList.addAll(1, this.removeSugarCandidates(tmpNewMolecule, tmpSugarCandidates));
        }
        //the molecule at index 0 may be empty and may be unconnected, based on the settings
        return tmpResultList;
    }




    public IAtomContainer removeCircularAndLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        IAtomContainer tmpNewMolecule = this.removeAndReturnCircularAndLinearSugars(aMolecule, aShouldBeCloned).get(0);
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    public boolean removeCircularAndLinearSugars(IAtomContainer aMolecule)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        List<IAtomContainer> tmpDeglycosylatedMoleculeAndSugarMoietiesList =
                this.removeAndReturnCircularAndLinearSugars(aMolecule, false);
        boolean tmpSomethingWasRemoved = (tmpDeglycosylatedMoleculeAndSugarMoietiesList.size() > 1);
        return tmpSomethingWasRemoved;
    }


    public List<IAtomContainer> removeAndReturnCircularAndLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        //throws NullPointerException if molecule is null
        this.setIndices(tmpNewMolecule);
        boolean tmpContainsCircularSugars = false;
        boolean tmpContainsLinearSugars = false;
        boolean tmpContainsAnyTypeOfSugars = false;
        //note: initial capacity arbitrarily chosen
        List<IAtomContainer> tmpResultList = new ArrayList<>(tmpNewMolecule.getAtomCount() / 6);
        tmpResultList.add(0, tmpNewMolecule);
        while (true) {
            //note: this has to be done stepwise because linear and circular sugar candidates can overlap
            //throws NullPointerException if molecule is null
            List<IAtomContainer> tmpCircularSugarCandidates = this.getCircularSugarCandidates(tmpNewMolecule);
            boolean tmpCandidateListIsNotEmpty = !tmpCircularSugarCandidates.isEmpty();
            List<IAtomContainer> tmpRemovedCircularSugarMoieties = new ArrayList<>(0);
            if (tmpCandidateListIsNotEmpty) {
                //throws NullPointerException and IllegalArgumentException
                tmpRemovedCircularSugarMoieties = this.removeSugarCandidates(tmpNewMolecule, tmpCircularSugarCandidates);
                if (!tmpContainsCircularSugars) {
                    tmpContainsCircularSugars = true;
                }
                tmpResultList.addAll(tmpRemovedCircularSugarMoieties);
            }
            //exit here if molecule is empty after removal
            if (tmpNewMolecule.isEmpty()) {
                break;
            }
            //note: if only terminal sugars are removed, the atom container should not be disconnected at this point
            // and that is a requirement for further checks for terminal linear sugar moieties
            //throws NullPointerException if molecule is null
            List<IAtomContainer> tmpLinearSugarCandidates = this.getLinearSugarCandidates(tmpNewMolecule);
            tmpCandidateListIsNotEmpty = !tmpLinearSugarCandidates.isEmpty();
            List<IAtomContainer> tmpRemovedLinearSugarMoieties = new ArrayList<>(0);
            if (tmpCandidateListIsNotEmpty) {
                //throws NullPointerException and IllegalArgumentException
                tmpRemovedLinearSugarMoieties = this.removeSugarCandidates(tmpNewMolecule, tmpLinearSugarCandidates);
                if (!tmpContainsLinearSugars) {
                    tmpContainsLinearSugars = true;
                }
                tmpResultList.addAll(tmpRemovedLinearSugarMoieties);
            }
            //exit here if molecule is empty after removal
            if (tmpNewMolecule.isEmpty()) {
                break;
            }
            if (this.removeOnlyTerminal) {
                int tmpCircularSugarCandidatesSizeAfterRemoval = tmpCircularSugarCandidates.size();
                int tmpLinearSugarCandidatesSizeAfterRemoval = tmpLinearSugarCandidates.size();
                boolean tmpSomethingWasRemoved = ((!tmpRemovedCircularSugarMoieties.isEmpty())
                        || (!tmpRemovedLinearSugarMoieties.isEmpty()));
                if (!tmpSomethingWasRemoved) {
                    //if nothing was removed, the loop is broken; otherwise, there might be new terminal moieties in the
                    // next iteration
                    break;
                }
            } else {
                break;
            }
        }
        if (this.setPropertyOfSugarContainingMolecules) {
            tmpContainsAnyTypeOfSugars = (tmpContainsCircularSugars || tmpContainsLinearSugars);
            tmpNewMolecule.setProperty(CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsAnyTypeOfSugars);
            tmpNewMolecule.setProperty(CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsCircularSugars);
            tmpNewMolecule.setProperty(CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsLinearSugars);
        }
        //The molecule at index 0 may be empty and may be unconnected, based on the settings
        return tmpResultList;
    }



    public static IAtomContainer selectBiggestUnconnectedFragment(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        if (tmpIsConnected) {
            return aMolecule;
        }
        Map<Object, Object> tmpProperties = aMolecule.getProperties();
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpBiggestFragment;
        if(tmpUnconnectedFragments != null && tmpUnconnectedFragments.getAtomContainerCount() >= 1) {
            tmpBiggestFragment = tmpUnconnectedFragments.getAtomContainer(0);
            int tmpBiggestFragmentHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(tmpBiggestFragment).size();
            for(IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()){
                int tmpFragmentHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(tmpFragment).size();
                if(tmpFragmentHeavyAtomCount > tmpBiggestFragmentHeavyAtomCount){
                    tmpBiggestFragment = tmpFragment;
                    tmpBiggestFragmentHeavyAtomCount = tmpFragmentHeavyAtomCount;
                }
            }
        } else {
            throw new NullPointerException("Could not detect the unconnected structures of the given atom container.");
        }
        tmpBiggestFragment.setProperties(tmpProperties);
        return tmpBiggestFragment;
    }


    public static IAtomContainer selectHeaviestUnconnectedFragment(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        if (tmpIsConnected) {
            return aMolecule;
        }
        Map<Object, Object> tmpProperties = aMolecule.getProperties();
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpHeaviestFragment;
        if(tmpUnconnectedFragments != null && tmpUnconnectedFragments.getAtomContainerCount() >= 1) {
            tmpHeaviestFragment = tmpUnconnectedFragments.getAtomContainer(0);
            double tmpHeaviestFragmentWeight = AtomContainerManipulator.getMass(tmpHeaviestFragment);
            for(IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()){
                double tmpFragmentWeight = AtomContainerManipulator.getMass(tmpFragment);
                if(tmpFragmentWeight > tmpHeaviestFragmentWeight){
                    tmpHeaviestFragment = tmpFragment;
                    tmpHeaviestFragmentWeight = tmpFragmentWeight;
                }
            }
        } else {
            //if something went wrong
            return null;
        }
        tmpHeaviestFragment.setProperties(tmpProperties);
        return tmpHeaviestFragment;
    }


    public static List<IAtomContainer> partitionAndSortUnconnectedFragments(IAtomContainer aMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        boolean tmpIsEmpty = aMolecule.isEmpty();
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        if (tmpIsConnected || tmpIsEmpty) {
            ArrayList<IAtomContainer> tmpFragmentList = new ArrayList<>(1);
            tmpFragmentList.add(aMolecule);
            return tmpFragmentList;
        }
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        int tmpSize = tmpUnconnectedFragments.getAtomContainerCount();
        ArrayList<IAtomContainer> tmpSortedList = new ArrayList<>(tmpSize);
        for (IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()) {
            tmpSortedList.add(tmpFragment);
        }
        /*Compares two IAtomContainers for order with the following criteria with decreasing priority:
            Compare atom count
            Compare molecular weight (heavy atoms only)
            Compare bond count
            Compare sum of bond orders (heavy atoms only)
        If no difference can be found with the above criteria, the IAtomContainers are considered equal.*/
        Comparator tmpComparator = new AtomContainerComparator().reversed();
        //note: this can throw various exceptions but they should not appear here
        tmpSortedList.sort(tmpComparator);
        return tmpSortedList;
    }





    protected void setIndices(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return;
        }
        for (int i = 0; i < aMolecule.getAtomCount(); i++) {
            IAtom tmpAtom = aMolecule.getAtom(i);
            tmpAtom.setProperty(this.INDEX_PROPERTY_KEY, i);
        }
    }


    protected Set<String> createSubstructureIdentifier(List<IAtomContainer> aSubstructureList) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aSubstructureList, "Given list is 'null'");
        for (IAtomContainer tmpSubstructure : aSubstructureList) {
            //method checks that no index appears multiple times but does not check whether there are numbers missing
            boolean tmpAreIndicesSet = this.checkUniqueIndicesOfAtoms(tmpSubstructure);
            if (!tmpAreIndicesSet) {
                throw new IllegalArgumentException("This method requires that every atom has a unique index.");
            }
        }
        HashSet<String> tmpIdentifierSet = new HashSet(aSubstructureList.size(), 1.0f);
        for (IAtomContainer tmpSubstructure: aSubstructureList) {
            if (Objects.isNull(tmpSubstructure)) {
                continue;
            }
            if (tmpSubstructure.isEmpty()) {
                continue;
            }
            String tmpSubstructureIdentifier = this.createSubstructureIdentifier(tmpSubstructure);
            tmpIdentifierSet.add(tmpSubstructureIdentifier);
        }
        return tmpIdentifierSet;
    }

    protected String createSubstructureIdentifier(IAtomContainer aSubstructure) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aSubstructure, "Given substructure is 'null'.");
        boolean tmpAreIndicesSet = this.checkUniqueIndicesOfAtoms(aSubstructure);
        if (!tmpAreIndicesSet) {
            throw new IllegalArgumentException("This method requires that every atom has a unique index.");
        }
        List<Integer> tmpIndicesList = new ArrayList<>(aSubstructure.getAtomCount());
        for (IAtom tmpAtom : aSubstructure.atoms()) {
            int tmpAtomIndex = tmpAtom.getProperty(this.INDEX_PROPERTY_KEY);
            tmpIndicesList.add(tmpAtomIndex);
        }
        Collections.sort(tmpIndicesList);
        String tmpSubstructureIdentifier = "";
        for (int tmpAtomIndex : tmpIndicesList) {
            tmpSubstructureIdentifier += Integer.toString(tmpAtomIndex);
        }
        return tmpSubstructureIdentifier;
    }


    private boolean areAllExocyclicBondsSingle(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule " +
                "is 'null'");
        int tmpAtomCountInRing = aRingToTest.getAtomCount();
        int tmpArrayListInitCapacity = tmpAtomCountInRing * 2;
        List<IBond> tmpExocyclicBondsList = new ArrayList<>(tmpArrayListInitCapacity);
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            List<IBond> tmpConnectedBondsList = anOriginalMolecule.getConnectedBondsList(tmpRingAtom);
            for (IBond tmpBond : tmpConnectedBondsList) {
                boolean tmpIsInRing = aRingToTest.contains(tmpBond);
                if (!tmpIsInRing) {
                    tmpExocyclicBondsList.add(tmpBond);
                }
            }
        }
        return (BondManipulator.getMaximumBondOrder(tmpExocyclicBondsList) == IBond.Order.SINGLE);
    }


    private boolean hasGlycosidicBond(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule " +
                "is 'null'");
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        boolean tmpContainsGlycosidicBond = false;
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            boolean tmpBreakOuterLoop = false;
            //check to avoid exceptions
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            List<IAtom> connectedAtomsList = anOriginalMolecule.getConnectedAtomsList(tmpRingAtom);
            for (IAtom tmpAtom : connectedAtomsList) {
                boolean tmpIsInRing = aRingToTest.contains(tmpAtom);
                if (!tmpIsInRing) {
                    String tmpSymbol = tmpAtom.getSymbol();
                    boolean tmpIsOxygen = (tmpSymbol.equals("O"));
                    if (tmpIsOxygen) {
                        List<IBond> tmpConnectedBondsList = anOriginalMolecule.getConnectedBondsList(tmpAtom);
                        boolean tmpHasOnlyTwoBonds = (tmpConnectedBondsList.size() == 2);
                        boolean tmpAllBondsAreSingle =
                                (BondManipulator.getMaximumBondOrder(tmpConnectedBondsList) == IBond.Order.SINGLE);
                        boolean tmpOneBondAtomIsHydrogen = false;
                        for (IBond tmpBond : tmpConnectedBondsList) {
                            for (IAtom tmpBondAtom : tmpBond.atoms()) {
                                if (tmpBondAtom.getSymbol().equals("H")) {
                                    tmpOneBondAtomIsHydrogen = true;
                                }
                            }
                        }
                        if ((tmpHasOnlyTwoBonds && tmpAllBondsAreSingle) && !tmpOneBondAtomIsHydrogen) {
                            tmpContainsGlycosidicBond = true;
                            tmpBreakOuterLoop = true;
                            break;
                        }
                    }
                }
            }
            if (tmpBreakOuterLoop) {
                break;
            }
        }
        return tmpContainsGlycosidicBond;
    }


    private int getAttachedOxygenAtomCount(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule " +
                "is 'null'");
        int tmpExocyclicOxygenCounter = 0;
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            //check to avoid exceptions
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            List<IAtom> tmpConnectedAtomsList = anOriginalMolecule.getConnectedAtomsList(tmpRingAtom);
            for (IAtom tmpConnectedAtom : tmpConnectedAtomsList) {
                String tmpSymbol = tmpConnectedAtom.getSymbol();
                boolean tmpIsOxygen = tmpSymbol.equals("O");
                boolean tmpIsInRing = aRingToTest.contains(tmpConnectedAtom);
                if (tmpIsOxygen && !tmpIsInRing) {
                    tmpExocyclicOxygenCounter++;
                }
            }
        }
        return tmpExocyclicOxygenCounter;
    }



    public void clearTooSmallStructures(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return;
        }
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            return;
        }
        IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        for (int i = 0; i < tmpComponents.getAtomContainerCount(); i++) {
            IAtomContainer tmpComponent = tmpComponents.getAtomContainer(i);
            //May throw UnsupportedOperationException if a new StructureToKeepMode option has been added but not implemented
            // in this method yet. Since this is a serious issue, the code is supposed to crash.
            boolean tmpIsTooSmall = this.isTooSmall(tmpComponent);
            if (tmpIsTooSmall) {
                //note: careful with removing things from sets/lists while iterating over it! But here it is ok because elements
                // are not removed from the same set that is iterated
                for (IAtom tmpAtom : tmpComponent.atoms()) {
                    //check to avoid exceptions
                    if (aMolecule.contains(tmpAtom)) {
                        aMolecule.removeAtom(tmpAtom);
                    }
                }
            }
        }
    }




    public boolean isTooSmall(IAtomContainer aMolecule) throws NullPointerException, UnsupportedOperationException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return true;
        }
        boolean tmpIsTooSmall;
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            tmpIsTooSmall = false;
        } else if (this.structuresToKeepMode == StructuresToKeepMode.HEAVY_ATOM_COUNT) {
            int tmpHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(aMolecule).size();
            tmpIsTooSmall = tmpHeavyAtomCount < this.structureToKeepModeThreshold;
        } else if (this.structuresToKeepMode == StructuresToKeepMode.MOLECULAR_WEIGHT) {
            double tmpMolWeight = AtomContainerManipulator.getMass(aMolecule, AtomContainerManipulator.MolWeight);
            tmpIsTooSmall = tmpMolWeight < this.structureToKeepModeThreshold;
        } else {
            throw new UnsupportedOperationException("Undefined StructuresToKeepMode setting!");
        }
        return tmpIsTooSmall;
    }


    public boolean isTerminal(IAtomContainer aSubstructure,
                              IAtomContainer aParentMolecule,
                              List<IAtomContainer> aCandidateList)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aSubstructure, "Given substructure is 'null'.");
        Objects.requireNonNull(aParentMolecule, "Given parent molecule is 'null'.");
        Objects.requireNonNull(aCandidateList, "Given list of candidates is 'null'.");
        boolean tmpIsParent = true;
        for (IAtom tmpAtom : aSubstructure.atoms()) {
            if (!aParentMolecule.contains(tmpAtom)) {
                tmpIsParent = false;
                break;
            }
        }
        if (!tmpIsParent) {
            throw new IllegalArgumentException("Given substructure is not part of the given parent molecule.");
        }
        boolean tmpIsUnconnected = !ConnectivityChecker.isConnected(aParentMolecule);
        if (tmpIsUnconnected) {
            throw new IllegalArgumentException("Parent molecule is already unconnected.");
        }
        boolean tmpIndicesAreSet = this.checkUniqueIndicesOfAtoms(aParentMolecule);
        if (!tmpIndicesAreSet) {
            this.setIndices(aParentMolecule);
        }
        //</editor-fold>
        boolean tmpIsTerminal;
        IAtomContainer tmpMoleculeClone = aParentMolecule.clone();
        IAtomContainer tmpSubstructureClone = aSubstructure.clone();
        HashMap<Integer, IAtom> tmpIndexToAtomMap = new HashMap<>(tmpMoleculeClone.getAtomCount() + 1, 1);
        for (IAtom tmpAtom : tmpMoleculeClone.atoms()) {
            tmpIndexToAtomMap.put(tmpAtom.getProperty(this.INDEX_PROPERTY_KEY), tmpAtom);
        }
        for (IAtom tmpAtom : tmpSubstructureClone.atoms()) {
            tmpMoleculeClone.removeAtom(tmpIndexToAtomMap.get(tmpAtom.getProperty(this.INDEX_PROPERTY_KEY)));
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(tmpMoleculeClone);
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            tmpIsTerminal = tmpIsConnected;
        } else {
            if (tmpIsConnected) {
                tmpIsTerminal = true;
            } else {
                IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(tmpMoleculeClone);
                HashSet<Integer> tmpAtomIndicesThatArePartOfSugarCandidates = new HashSet<>(aParentMolecule.getAtomCount(), 0.8f);
                for (IAtomContainer tmpCandidate : aCandidateList) {
                    for (IAtom tmpAtom : tmpCandidate.atoms()) {
                        int tmpIndex = tmpAtom.getProperty(this.INDEX_PROPERTY_KEY);
                        tmpAtomIndicesThatArePartOfSugarCandidates.add(tmpIndex);
                    }
                }
                for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
                    if (Objects.isNull(tmpComponent) || tmpComponent.isEmpty()) {
                        continue;
                    }
                    //May throw UnsupportedOperationException if a new StructureToKeepMode option has been added but not implemented
                    // in this method yet. Since this is a serious issue, the code is supposed to crash.
                    //throws NullPointerException if molecule is null
                    boolean tmpIsTooSmall = this.isTooSmall(tmpComponent);
                    boolean tmpIsPartOfSugarCandidate = false;
                    for (IAtom tmpAtom : tmpComponent.atoms()) {
                        int tmpIndex = tmpAtom.getProperty(this.INDEX_PROPERTY_KEY);
                        if (tmpAtomIndicesThatArePartOfSugarCandidates.contains(tmpIndex)) {
                            tmpIsPartOfSugarCandidate = true;
                            break;
                        }
                    }
                    if (tmpIsTooSmall && !tmpIsPartOfSugarCandidate) {
                        //note: no check whether the clone actually contains the component
                        tmpMoleculeClone.remove(tmpComponent);
                    }
                }
                tmpIsTerminal = ConnectivityChecker.isConnected(tmpMoleculeClone);
            }
        }
        return tmpIsTerminal;
    }

    protected boolean checkUniqueIndicesOfAtoms(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Given molecule is empty.");
        }
        HashSet<Integer> tmpAtomIndices = new HashSet<>(aMolecule.getAtomCount() + 4, 1.0f);
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (Objects.isNull(tmpAtom.getProperty(this.INDEX_PROPERTY_KEY))) {
                return false;
            } else {
                int tmpIndex = tmpAtom.getProperty(this.INDEX_PROPERTY_KEY);
                if (tmpAtomIndices.contains(tmpIndex)) {
                    return false;
                } else {
                    tmpAtomIndices.add(tmpIndex);
                }
            }
        }
        //only reached if method is not exited before because of a missing or non-unique index
        return true;
    }



    public void postProcessAfterRemoval(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        if (aMolecule.isEmpty()) {
            return;
        }
        //if too small / too light, unconnected structures should be discarded, this is done now
        //otherwise, the possibly unconnected atom container is returned
        //Even if only terminal sugars are removed, the resulting, connected structure may still be too small to keep!
        if (this.structuresToKeepMode != StructuresToKeepMode.ALL) {
            //throws NullPointerException if molecule is null
            this.clearTooSmallStructures(aMolecule);
        }
        if (!aMolecule.isEmpty()) {
            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
                CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(aMolecule);
            } catch (CDKException aCDKException) {
                aCDKException.printStackTrace();
            }
        }
    }



    public List<IAtomContainer> getCircularSugarCandidates(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        if (aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        boolean tmpIndicesAreSet = this.checkUniqueIndicesOfAtoms(aMolecule);
        if (!tmpIndicesAreSet) {
            this.setIndices(aMolecule);
        }
        int[][] tmpAdjList = GraphUtil.toAdjList(aMolecule);
        //efficient computation/partitioning of the ring systems
        RingSearch tmpRingSearch = new RingSearch(aMolecule, tmpAdjList);
        List<IAtomContainer> tmpPotentialSugarRings = this.getPotentialSugarCycles(aMolecule);
        if (tmpPotentialSugarRings.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(tmpPotentialSugarRings.size());
        for(IAtomContainer tmpPotentialSugarRing : tmpPotentialSugarRings) {
            if (Objects.isNull(tmpPotentialSugarRing) || tmpPotentialSugarRing.isEmpty()) {
                continue;
            }
            /*note: another requirement of a suspected sugar ring is that it contains only single bonds.
             * This is not tested here because all the structures in the reference rings do meet this criterion.
             * But a structure that does not meet this criterion could be added to the references by the user.*/
            //do not remove rings without an attached glycosidic bond if this option is set
            if (this.detectGlycosidicBond) {
                boolean tmpHasGlycosidicBond = this.hasGlycosidicBond(tmpPotentialSugarRing, aMolecule);
                if (!tmpHasGlycosidicBond) {
                    //special exemption for molecules that only consist of a sugar ring and nothing else:
                    // they should also be seen as candidate even though they do not have a glycosidic bond
                    // (because there is nothing to bind to)
                    if (tmpRingSearch.numRings() == 1) {
                        boolean tmpMoleculeIsOnlyOneSugarRing = false;
                        try {
                            tmpMoleculeIsOnlyOneSugarRing = this.checkCircularSugarGlycosidicBondExemption(tmpPotentialSugarRing, aMolecule);
                        } catch (CloneNotSupportedException | IllegalArgumentException | NullPointerException anException) {
                            anException.printStackTrace();
                            //there is sth wrong here, do not add this ring to the candidates
                            continue;
                        }
                        if (!tmpMoleculeIsOnlyOneSugarRing) {
                            //isolated ring is not a candidate because it has no glycosidic bond and does not
                            // qualify for the exemption
                            continue;
                        } //else, go on investigating this candidate, even though it does not have a glycosidic bond
                    } else {
                        //not a candidate
                        continue;
                    }
                }
            }
            //do not remove rings with 'too few' attached oxygens if this option is set
            if (this.includeNrOfAttachedOxygens) {
                int tmpExocyclicOxygenCount = this.getAttachedOxygenAtomCount(tmpPotentialSugarRing, aMolecule);
                int tmpAtomsInRing = tmpPotentialSugarRing.getAtomCount();
                boolean tmpAreEnoughOxygensAttached = this.doesRingHaveEnoughOxygenAtomsAttached(tmpAtomsInRing,
                        tmpExocyclicOxygenCount);
                if (!tmpAreEnoughOxygensAttached) {
                    continue;
                }
            }
            //if sugar ring has not been excluded yet, the molecule contains sugars, although they might not
            // be terminal
            tmpSugarCandidates.add(tmpPotentialSugarRing);
        }
        return tmpSugarCandidates;
    }

    protected List<IAtomContainer> getPotentialSugarCycles(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        boolean tmpIndicesAreSet = this.checkUniqueIndicesOfAtoms(aMolecule);
        if (!tmpIndicesAreSet) {
            this.setIndices(aMolecule);
        }
        int[][] tmpAdjList = GraphUtil.toAdjList(aMolecule);
        //efficient computation/partitioning of the ring systems
        RingSearch tmpRingSearch = new RingSearch(aMolecule, tmpAdjList);
        List<IAtomContainer> tmpIsolatedRings = tmpRingSearch.isolatedRingFragments();
        if (tmpIsolatedRings.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(tmpIsolatedRings.size());
        for(IAtomContainer tmpReferenceRing : this.ringSugars) {
            for (IAtomContainer tmpIsolatedRing : tmpIsolatedRings) {
                if (Objects.isNull(tmpIsolatedRing) || tmpIsolatedRing.isEmpty()) {
                    continue;
                }
                boolean tmpIsIsomorph = false;
                UniversalIsomorphismTester tmpUnivIsoTester = new UniversalIsomorphismTester();
                try {
                    tmpIsIsomorph = tmpUnivIsoTester.isIsomorph(tmpReferenceRing, tmpIsolatedRing);
                } catch (CDKException aCDKException) {
                    aCDKException.printStackTrace();
                    continue;
                }
                if (tmpIsIsomorph) {
                    /*note: another requirement of a suspected sugar ring is that it contains only single bonds.
                     * This is not tested here because all the structures in the reference rings do meet this criterion.
                     * But a structure that does not meet this criterion could be added to the references by the user.*/
                    //do not remove rings with non-single exocyclic bonds, they are not sugars (not an option!)
                    boolean tmpAreAllExocyclicBondsSingle = this.areAllExocyclicBondsSingle(tmpIsolatedRing, aMolecule);
                    if (tmpAreAllExocyclicBondsSingle) {
                        tmpSugarCandidates.add(tmpIsolatedRing);
                    }
                }
            }
        }
        return tmpSugarCandidates;
    }

    protected boolean checkCircularSugarGlycosidicBondExemption(IAtomContainer aRing, IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aRing, "Given ring is 'null'.");
        Objects.requireNonNull(aMolecule, "Given parent molecule is 'null'.");
        boolean tmpIsParent = true;
        for (IAtom tmpAtom : aRing.atoms()) {
            if (!aMolecule.contains(tmpAtom)) {
                tmpIsParent = false;
                break;
            }
        }
        if (!tmpIsParent) {
            throw new IllegalArgumentException("Given substructure is not part of the given parent molecule.");
        }
        boolean tmpIndicesAreSet = this.checkUniqueIndicesOfAtoms(aMolecule);
        if (!tmpIndicesAreSet) {
            this.setIndices(aMolecule);
        }
        boolean tmpQualifiesForExemption = false;
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        IAtomContainer tmpSubstructureClone = aRing.clone();
        HashMap<Integer, IAtom> tmpIndexToAtomMap = new HashMap<>(tmpMoleculeClone.getAtomCount() + 1, 1.0f);
        for (IAtom tmpAtom : tmpMoleculeClone.atoms()) {
            tmpIndexToAtomMap.put(tmpAtom.getProperty(INDEX_PROPERTY_KEY), tmpAtom);
        }
        for (IAtom tmpAtom : tmpSubstructureClone.atoms()) {
            tmpMoleculeClone.removeAtom(tmpIndexToAtomMap.get(tmpAtom.getProperty(INDEX_PROPERTY_KEY)));
        }
        if (tmpMoleculeClone.isEmpty()) {
            tmpQualifiesForExemption = true;
        } else {
            this.clearTooSmallStructures(tmpMoleculeClone);
            tmpQualifiesForExemption = tmpMoleculeClone.isEmpty();
        }
        return tmpQualifiesForExemption;
    }




    public List<IAtomContainer> getLinearSugarCandidates(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        if (aMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        boolean tmpIndicesAreSet = this.checkUniqueIndicesOfAtoms(aMolecule);
        if (!tmpIndicesAreSet) {
            this.setIndices(aMolecule);
        }
        List<IAtomContainer> tmpSugarCandidates = this.linearSugarCandidatesByPatternMatching(aMolecule);
        //alternative: SMARTS or Ertl or matching the biggest patterns first and exclude the matched atoms
        if (!tmpSugarCandidates.isEmpty()) {

            //*Debugging*
            //this.printAllMolsAsSmiles(tmpSugarCandidates);

            tmpSugarCandidates = this.combineOverlappingCandidates(tmpSugarCandidates);
            //alternative: tmpSugarCandidates = this.splitOverlappingCandidatesPseudoRandomly(tmpSugarCandidates);

            //*Debugging*
            //this.printAllMolsAsSmiles(tmpSugarCandidates);

            tmpSugarCandidates = this.splitEtherEsterAndPeroxideBonds(tmpSugarCandidates);

            //*Debugging*
            //this.printAllMolsAsSmiles(tmpSugarCandidates);

            this.removeCandidatesContainingCircularSugars(tmpSugarCandidates, aMolecule);
            //alternative: tmpSugarCandidates = this.removeCircularSugarsFromCandidates(tmpSugarCandidates);

            //*Debugging*
            //this.printAllMolsAsSmiles(tmpSugarCandidates);

            tmpSugarCandidates = this.removeTooSmallAndTooLargeCandidates(tmpSugarCandidates);

            //*Debugging*
            //this.printAllMolsAsSmiles(tmpSugarCandidates);
        }
        if (!this.removeLinearSugarsInRing && !tmpSugarCandidates.isEmpty()) {
            this.removeSugarCandidatesWithCyclicAtoms(tmpSugarCandidates, aMolecule);
            //alternative: tmpSugarCandidates = this.removeCyclicAtomsFromSugarCandidates(tmpSugarCandidates, tmpNewMolecule);
        }
        return tmpSugarCandidates;
    }



    protected void removeSugarCandidatesWithCyclicAtoms(List<IAtomContainer> aCandidateList,
                                                        IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty()) {
            return;
        }
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int[][] tmpAdjList = GraphUtil.toAdjList(aMolecule);
        RingSearch tmpRingSearch = new RingSearch(aMolecule, tmpAdjList);
        for (int i = 0; i < aCandidateList.size(); i++) {
            IAtomContainer tmpCandidate = aCandidateList.get(i);
            for (int j = 0; j < tmpCandidate.getAtomCount(); j++) {
                IAtom tmpAtom = tmpCandidate.getAtom(j);
                if (tmpRingSearch.cyclic(tmpAtom)) {
                    aCandidateList.remove(i);
                    //removal shifts the remaining indices
                    i = i - 1;
                    break;
                }
            }
        }
    }


    protected List<IAtomContainer> removeTooSmallAndTooLargeCandidates(List<IAtomContainer> aCandidateList) throws NullPointerException {
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty()) {
            return aCandidateList;
        }
        List<IAtomContainer> tmpProcessedCandidates = new ArrayList<>(aCandidateList.size());
        for (IAtomContainer tmpCandidate : aCandidateList) {
            int tmpCarbonCount = 0;
            for (IAtom tmpAtom : tmpCandidate.atoms()) {
                String tmpSymbol = tmpAtom.getSymbol();
                if (tmpSymbol.equals("C")) {
                    tmpCarbonCount++;
                }
            }
            if (tmpCarbonCount >= this.linearSugarCandidateMinSize && tmpCarbonCount <= this.linearSugarCandidateMaxSize) {
                tmpProcessedCandidates.add(tmpCandidate);
            }
        }
        return tmpProcessedCandidates;
    }



    protected void removeCandidatesContainingCircularSugars(List<IAtomContainer> aCandidateList,
                                                            IAtomContainer aParentMolecule)
            throws NullPointerException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty()) {
            return;
        }
        Objects.requireNonNull(aParentMolecule, "Given parent molecule is 'null'.");
        //</editor-fold>
        // generating ids for the isolated potential sugar circles in the parent molecule
        List<IAtomContainer> tmpPotentialSugarRingsParent = this.getPotentialSugarCycles(aParentMolecule);
        // nothing to process
        if (tmpPotentialSugarRingsParent.isEmpty()) {
            return;
        }
        Set<String> tmpPotentialSugarRingsParentIdentifierSet = this.createSubstructureIdentifier(tmpPotentialSugarRingsParent);
        // iterating over candidates
        for (int i = 0; i < aCandidateList.size(); i++) {
            IAtomContainer tmpCandidate = aCandidateList.get(i);
            if (Objects.isNull(tmpCandidate)) {
                aCandidateList.remove(i);
                //The removal shifts the remaining indices!
                i = i - 1;
                continue;
            }
            List<IAtomContainer> tmpPotentialSugarRingsCandidate = this.getPotentialSugarCycles(tmpCandidate);
            boolean tmpIsAlsoIsolatedInParent = false;
            if (!tmpPotentialSugarRingsCandidate.isEmpty()) {
                //iterating over potential sugar rings in candidate
                for(IAtomContainer tmpRing : tmpPotentialSugarRingsCandidate) {
                    if (Objects.isNull(tmpRing) || tmpRing.isEmpty()) {
                        continue;
                    }
                    String tmpRingIdentifier = this.createSubstructureIdentifier(tmpRing);
                    tmpIsAlsoIsolatedInParent = tmpPotentialSugarRingsParentIdentifierSet.contains(tmpRingIdentifier);
                    if (tmpIsAlsoIsolatedInParent) {
                        aCandidateList.remove(i);
                        i = i - 1;
                        // break the iteration of rings and go to the next candidate
                        break;
                    }
                }
            }
        }
    }



    protected List<IAtomContainer> splitEtherEsterAndPeroxideBonds(List<IAtomContainer> aCandidateList) throws NullPointerException {
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        List<IAtomContainer> tmpProcessedCandidates = new ArrayList<>(aCandidateList.size() * 2);
        for (IAtomContainer tmpCandidate : aCandidateList) {
            SmartsPattern.prepare(tmpCandidate);

            // note: ester matching has to precede the ether matching because the ether pattern also matches esters
            // note 2: here, which bond is removed is specifically defined. This is not the case for the ether
            Mappings tmpEsterMappings = this.ESTER_SMARTS_PATTERN.matchAll(tmpCandidate).uniqueAtoms();
            if (tmpEsterMappings.atLeast(1)) {
                for (IAtomContainer tmpEsterGroup : tmpEsterMappings.toSubstructures()) {
                    IAtom tmpDoubleBondedOxygen = null;
                    IAtom tmpConnectingOxygen = null;
                    for (IAtom tmpAtom : tmpEsterGroup.atoms()) {
                        String tmpSymbol = tmpAtom.getSymbol();
                        if (tmpSymbol.equals("O")) {
                            int tmpBondCount = tmpAtom.getBondCount();
                            if (tmpBondCount == 1) {
                                tmpDoubleBondedOxygen = tmpAtom;
                            } else {
                                tmpConnectingOxygen = tmpAtom;
                            }
                        }
                    }
                    IAtom tmpCarbonBoundToDoubleBondedOxygen = tmpEsterGroup.getConnectedAtomsList(tmpDoubleBondedOxygen).get(0);
                    tmpCandidate.removeBond(tmpCarbonBoundToDoubleBondedOxygen, tmpConnectingOxygen);
                }
            }

            // note: which bond is actually removed is 'random', i.e. not to predict by a human
            Mappings tmpEtherMappings = this.ETHER_SMARTS_PATTERN.matchAll(tmpCandidate).uniqueAtoms();
            if (tmpEtherMappings.atLeast(1)) {
                for (IAtomContainer tmpEtherGroup : tmpEtherMappings.toSubstructures()) {
                    IAtom tmpCarbon1 = null;
                    IAtom tmpCarbon2 = null;
                    IAtom tmpOxygen = null;
                    for (IAtom tmpAtom : tmpEtherGroup.atoms()) {
                        String tmpSymbol = tmpAtom.getSymbol();
                        if (tmpSymbol.equals("O")) {
                            tmpOxygen = tmpAtom;
                        } else if (tmpSymbol.equals("C") && Objects.isNull(tmpCarbon1)) {
                            tmpCarbon1 = tmpAtom;
                        } else {
                            tmpCarbon2 = tmpAtom;
                        }
                    }
                    tmpCandidate.removeBond(tmpOxygen, tmpCarbon2);
                }
            }

            Mappings tmpPeroxideMappings = this.PEROXIDE_SMARTS_PATTERN.matchAll(tmpCandidate).uniqueAtoms();
            if (tmpPeroxideMappings.atLeast(1)) {
                for (IAtomContainer tmpPeroxideGroup : tmpPeroxideMappings.toSubstructures()) {
                    IAtom tmpOxygen1 = null;
                    IAtom tmpOxygen2 =  null;
                    for (IAtom tmpAtom : tmpPeroxideGroup.atoms()) {
                        String tmpSymbol = tmpAtom.getSymbol();
                        if (tmpSymbol.equals("O")) {
                            if (Objects.isNull(tmpOxygen1)) {
                                tmpOxygen1 = tmpAtom;
                            } else {
                                tmpOxygen2 = tmpAtom;
                            }
                        }
                    }
                    tmpCandidate.removeBond(tmpOxygen1, tmpOxygen2);
                }
            }

            boolean tmpIsConnected = ConnectivityChecker.isConnected(tmpCandidate);
            if (tmpIsConnected) {
                tmpProcessedCandidates.add(tmpCandidate);
            } else {
                IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(tmpCandidate);
                for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
                    tmpProcessedCandidates.add(tmpComponent);
                }
            }
        }
        return tmpProcessedCandidates;
    }



    protected List<IAtomContainer> linearSugarCandidatesByPatternMatching(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        IAtomContainer tmpNewMolecule = aMolecule;
        if (tmpNewMolecule.isEmpty()) {
            return new ArrayList<IAtomContainer>(0);
        }
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(tmpNewMolecule.getAtomCount() / 2);
        for (DfPattern tmpLinearSugarPattern : this.linearSugarPatterns) {
            if (Objects.isNull(tmpLinearSugarPattern)) {
                continue;
            }
            /*unique in this case means that the same match cannot be in this collection multiple times but they can
            still overlap!*/
            Mappings tmpMappings = tmpLinearSugarPattern.matchAll(tmpNewMolecule);
            Mappings tmpUniqueMappings = tmpMappings.uniqueAtoms();
            Iterable<IAtomContainer> tmpUniqueSubstructureMappings = tmpUniqueMappings.toSubstructures();
            for (IAtomContainer tmpMatchedStructure : tmpUniqueSubstructureMappings) {
                if (Objects.isNull(tmpMatchedStructure)) {
                    continue;
                }
                tmpSugarCandidates.add(tmpMatchedStructure);
            }
        }
        return tmpSugarCandidates;
    }


    protected List<IAtomContainer> combineOverlappingCandidates(List<IAtomContainer> aCandidateList) throws NullPointerException  {
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty()) {
            return aCandidateList;
        }
        int tmpListSize = aCandidateList.size();
        List<IAtomContainer> tmpNonOverlappingSugarCandidates = new ArrayList<>(tmpListSize);
        IAtomContainer tmpMatchesContainer = new AtomContainer();
        for (int i = 0; i < tmpListSize; i++) {
            IAtomContainer tmpCandidate = aCandidateList.get(i);
            tmpMatchesContainer.add(tmpCandidate);
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(tmpMatchesContainer);
        if (tmpIsConnected) {
            tmpNonOverlappingSugarCandidates.add(tmpMatchesContainer);
        } else {
            IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(tmpMatchesContainer);
            Iterable<IAtomContainer> tmpMolecules = tmpComponents.atomContainers();
            for (IAtomContainer tmpComponent : tmpMolecules) {
                tmpNonOverlappingSugarCandidates.add(tmpComponent);
            }
        }
        return tmpNonOverlappingSugarCandidates;
    }






    private boolean doesRingHaveEnoughOxygenAtomsAttached(int aNumberOfAtomsInRing,
                                                          int aNumberOfAttachedExocyclicOxygenAtoms) {
        if (aNumberOfAtomsInRing == 0) {
            //better than throwing an exception here?
            return false;
        }
        double tmpAttachedOxygensToAtomsInRingRatio =
                ((double) aNumberOfAttachedExocyclicOxygenAtoms / (double) aNumberOfAtomsInRing);
        boolean tmpMeetsThreshold =
                (tmpAttachedOxygensToAtomsInRingRatio >= this.attachedOxygensToAtomsInRingRatioThreshold);
        return tmpMeetsThreshold;
    }



    public void setAttachedOxygensToAtomsInRingRatioThreshold(double aDouble) throws IllegalArgumentException {
        //false for NaN and infinity arguments
        boolean tmpIsFinite = Double.isFinite(aDouble);
        boolean tmpIsNegative = (aDouble < 0);
        if(!tmpIsFinite || tmpIsNegative) {
            throw new IllegalArgumentException("Given double is NaN, infinite or negative.");
        }
        if (!this.includeNrOfAttachedOxygens) {
            throw new IllegalArgumentException("The number of attached oxygen atoms is currently not included in the " +
                    "decision making process, so a ratio threshold makes no sense.");
        }
        this.attachedOxygensToAtomsInRingRatioThreshold = aDouble;
    }


    public void setStructuresToKeepMode(StructuresToKeepMode aMode) throws NullPointerException {
        Objects.requireNonNull(aMode, "Given mode is 'null'.");
        this.structuresToKeepMode = aMode;
        this.structureToKeepModeThreshold = this.structuresToKeepMode.getDefaultThreshold();
    }

    public void setStructuresToKeepThreshold(int aThreshold) throws IllegalArgumentException {

        if ((this.structuresToKeepMode == StructuresToKeepMode.ALL)) {
            throw new IllegalArgumentException("The mode is currently set to keep all structures, so a threshold " +
                    "makes no sense.");
        }
        if (aThreshold < 0) {
            throw new IllegalArgumentException("Threshold cannot be negative.");
        }
        this.structureToKeepModeThreshold = aThreshold;
    }


}
