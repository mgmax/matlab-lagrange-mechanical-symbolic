Obtain symbolic ODEs for a mechanical system using the Lagrange method.
The code is from my bachelor's thesis, so use it with some caution, or rather as an inspiration to recreate it better.

For MATLAB R2013a (or newer, must have MuPad).

The remaining documentation is in German (sorry).


Mechanische Modellbildung für eine Dreh-Schwenk-Einheit:

mechanik_symbolisch.m öffnen und starten. Es erscheint ein MuPAD-Fenster, dort kann man mit den Gleichungen spielen, siehe dort enthaltene Beispiele

ACHTUNG: Beim ersten Start des Mupad-Fensters erscheinen die Gleichungen noch unschön (viele Variablen zu sigma_n zusammengefasst).
Deshalb in MuPad : View->Configure->Format->Output Settings->Abbreviate Output AUS schalten, Fenster schließen (nicht speichern) und Skript nochmal ausführen

Außerdem werden auch die tex-Ausgaben erzeugt (*.txt im dortigen Verzeichnis), die Anweisungen dazu stehen gegen Ende des Skriptes.

