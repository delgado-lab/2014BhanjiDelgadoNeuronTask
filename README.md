**Description of Task Files for "PersTaskfmri"**

**Used in publications**:

> Bhanji & Delgado (2014), also see Bhanji, Kim, & Delgado (2016), Bhanji, Delgado, & Ray (2021)

**Important Filenames (in the TaskFiles folder):**

> *PersTask_fmri_InstructPractice* - e-prime, instructions and practice
> before entering scanner
>
> *PersTaskfmri_order1* - e-prime, randomized order used for half of
> participants
>
> *PersTaskfmri_order2* - e-prime, randomized order used for half of
> participants- unc/con blocks switched from order1
>
> *images* - folder of images needed for the e-prime presentation
>
> *read_persist_logfiles.m* -- Matlab script, reads e-prime output files
> and generates FSL style 3-column timing files

**File Descriptions**

***PersTask_fmri_InstructPractice:***

> Instructions for task framed as "Academic degree game", followed by 2
> rounds of practice where the goal is not reached

***PersTaskfmri_order1 and PersTaskfmri_order2:***

*Experimental Design:* 2x2 within subjects (Condition changes from round
to round, but not within)

> Factor 1: Uncontrollable Obstacles (setbacks framed as random) versus
> Controllable Obstacles (setbacks framed as due to incorrect response)
> -- participant receives same pattern of setbacks in both conditions.
>
> Factor 2: High Alternative Value (Path Values are 80/78/76) versus Low
> Alternative Value (Path Values are 80/70/60)

*Timing Information:*

> Structured for 4 equal length (10m 30s) scanning runs (break screen
> between each)
>
> Event Timing - 2s Path Choice, 2/4/6s (50/25/25%) Fixation, 2s
> Obstacle Cue (includes response), 2/4/6s (50/25/25%) Fixation, 2s
> Obstacle Outcome (Setback received/avoided), 2/4/6s (50/25/25%)
> Fixation (\*no fixation between cue and outcome for Progress Cues)

*Trial Counts and other details:*

> 40 rounds, 128 obstacle cues, 80 setbacks (62.5% of obstacles), 48
> avoided setbacks (37.5% of obstacles), 64 progress cues (class
> meetings) -- divided equally across 4 conditions

*How to interpret output files:*

> EventType: 1=uncontrollable obstacle, 2=controllable obstacle,
> 3=progress cue, 4=path choice, 6=goal feedback
>
> Lose: 1=setback received, 0=setback avoided, -1=event was path choice
> or goal feedback
>
> Persist: 1 = choice to try again on the same path where a setback was
> just experienced, 0 = choice to switch to a different path than the
> one where the setback was experienced, -2=first path choice of the
> round, not included in persistence calculation, -1 = no response given
> for path choice

**read_persist_logfiles.m:** See documentation in script.

**Also see (e.g., similar tasks described elsewhere):** Similar to
PersTask_UncCon (used in StressPersist), PathTask (used in SmokPersist
and Opiod user studies)

**Other notes (e.g., how to calculate behavioral measures, other
versions of the task that might be helpful)**

> To calculate Behavioral Persistence for each participant, you take the
> number of Persist choices (coded as "1" in the "Persist" column of
> output for each condition, and divide by the total number of
> post-setback choices (choices with a "1" or "0" in the Persist column.
> Missed (no reponse) choices have a "-2" in the "Persist" column).
>
> The "Switch" column in the output has value "1" when the participant
> chooses a path that is lower in value than the current path, "0" if
> they choose the same path or a higher value path.
>
> If a participant does not change their response after a controllable
> setback (i.e., fails to learn from a mistake), then a setback is
> received and the game continues on.
>
> If a participant does not choose a Path in the given time then the
> highest value path is selected (for first choice in round) or their
> previous choice is selected (for subsequent choices), and the study
> continues on (these choices are marked with a value of "-2" in the
> "Persist" output column (indicating a missed response).
> 
> Other versions of the task are located in repositories in this 
> profile - for [Bhanji, Delgado, & Ray (2021](https://github.com/delgado-lab/POUD-PersistenceTask)
> and [Bhanji, Kim, Delgado (2016)(https://github.com/delgado-lab/PersistTask-JEPG).
>

