from PyQt5.QtWidgets import QWidget, QHBoxLayout, QLabel, QApplication, QVBoxLayout, QPushButton
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from io import BytesIO

def show_molecule_comparison(smiles1, smiles2):
    """
    Show a window comparing two molecules side by side.

    Args:
        smiles1: SMILES string for first molecule
        smiles2: SMILES string for second molecule

    Returns:
        tuple: (app, widget) - Keep reference to app if needed
    """
    # Create QApplication if it doesn't exist
    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    widget = MoleculeComparisonWidget(smiles1, smiles2)
    widget.setWindowTitle("CobraMod")
    widget.show()

    return app, widget

class MoleculeComparisonWidget(QWidget):
    def __init__(self, smiles1, smiles2, size=(400, 400)):
        super().__init__()
        self.size = size

        # Main layout
        self.main_layout = QVBoxLayout(self)

        # Info tables layout
        info_layout = QHBoxLayout()

        # Left info table
        left_info = QLabel()
        left_info.setText(
            "<b>Molecule 1</b><br>"
            "<table cellpadding='3'>"
            "<tr><td><b>Original Database:</b></td><td>KEGG</td></tr>"
            "<tr><td><b>Database ID:</b></td><td>C00031</td></tr>"
            "<tr><td><b>SMILES:</b></td><td>CC(=O)C(=O)O</td></tr>"
            "<tr><td><b>InChI:</b></td><td>InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)</td></tr>"
            "<tr><td><b>InChIKey:</b></td><td>LCTONWCANYUPML-UHFFFAOYSA-N</td></tr>"
            "<tr><td><b>Chemical Formular:</b></td><td>td></tr>"
            "</table>"
        )
        left_info.setWordWrap(True)
        left_info.setAlignment(Qt.AlignmentFlag.AlignTop)

        # Right info table
        right_info = QLabel()
        right_info.setText(
            "<b>Molecule 2</b><br>"
            "<table cellpadding='3'>"
            "<tr><td><b>Original Database:</b></td><td>BioCyc</td></tr>"
            "<tr><td><b>Database ID:</b></td><td>PYRUVATE</td></tr>"
            "<tr><td><b>SMILES:</b></td><td>CC(=O)C(=O)O</td></tr>"
            "<tr><td><b>InChI:</b></td><td>InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)</td></tr>"
            "<tr><td><b>InChIKey:</b></td><td>LCTONWCANYUPML-UHFFFAOYSA-N</td></tr>"
            "<tr><td><b>Chemical Formular:</b></td><td></td></tr>"
            "</table>"
        )
        right_info.setWordWrap(True)
        right_info.setAlignment(Qt.AlignmentFlag.AlignTop)

        info_layout.addWidget(left_info)
        info_layout.addWidget(right_info)

        self.main_layout.addLayout(info_layout)

        # Molecule display layout
        mol_layout = QHBoxLayout()

        # Create two labels for molecules
        self.label1 = QLabel()
        self.label2 = QLabel()

        self.label1.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label2.setAlignment(Qt.AlignmentFlag.AlignCenter)

        mol_layout.addWidget(self.label1)
        mol_layout.addWidget(self.label2)


        # Button layout
        button_layout = QHBoxLayout()

        btn_identical = QPushButton("Identical (Add once)")
        btn_different = QPushButton("Non Identical (Add Both)")

        btn_identical.clicked.connect(lambda: self.on_button_click("identical"))
        btn_different.clicked.connect(lambda: self.on_button_click("different"))

        button_layout.addWidget(btn_identical)
        button_layout.addWidget(btn_different)

        self.main_layout.addLayout(mol_layout)
        self.main_layout.addLayout(button_layout)

        # Display molecules
        self.set_molecules(smiles1, smiles2)

    def smiles_to_pixmap(self, smiles):
        """Convert SMILES to QPixmap"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=self.size)

        buffer = BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)

        qimg = QImage.fromData(buffer.getvalue())
        return QPixmap.fromImage(qimg)

    def set_molecules(self, smiles1, smiles2):
        """Update the displayed molecules"""
        pixmap1 = self.smiles_to_pixmap(smiles1)
        pixmap2 = self.smiles_to_pixmap(smiles2)

        if pixmap1:
            self.label1.setPixmap(pixmap1)
        else:
            self.label1.setText("Invalid SMILES 1")

        if pixmap2:
            self.label2.setPixmap(pixmap2)
        else:
            self.label2.setText("Invalid SMILES 2")
