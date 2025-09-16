# data/documentation/orig_names/yws_clueexcel/

This subdirectory contains the code & original output csvs from converting the Excel sheet with 12 different dataframes into 12 excel files (technically 13, since I did not exclude "Cover" in the python script. However, I did not copy that csv into `./data/yws_clue/`, since it is just documentation and would be nonsensical.)*

* Python script was run in the venv (Recall: source venv/bin/activate; deactivate; which python; which python3; which pip; *run these commands in Terminal to make sure I am installing packages & running the python script in a virtual environment, straight from Terminal (e.g. python3 yws_cluexcel.py)*)

* KEY: Translation of output names in this folder (which correpond directly to the sheet names in the Excel) to the csv names used in the code for the `menst_mod` repo as a whole, i.e. the names in `./data/yws_clue`: "output_..."

	- ALCOHOL_TOBACCO.csv = `ywsclue_alctob.csv`
	- CLUE_Birth_Control.csv = `ywsclue_bc.csv`
	- CLUE_Participants.csv = `ywsclue_part.csv`
	- CLUE_Sessions.csv = `ywsclue_sess.csv`
	- CLUE_Tracking.csv = `ywsclue_track.csv`

	- Derived.csv = `ywsclue_der.csv`
	- MEDICAL_HISTORY.csv = `ywsclue_medhist.csv`
	- Occupational_History.csv = `ywsclue_occhist.csv`
	- PHYSICAL_ACTIVITY.csv = `ywsclue_pa.csv`
	- REPRODUCTIVE.csv = `ywsclue_repro.csv`

	- Sleep.csv = `ywsclue_sleep.csv`
	- STRAIN.csv = `ywsclue_strain.csv`

