package ge.dchechelashvili.notesapp.main

import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

interface IMainView {
    abstract fun showPinnedNotes(pinnedNotes: List<Note>, pinnedItems: List<List<Item>>)
    abstract fun showUnpinnedNotes(unpinnedNotes: List<Note>, unpinnedItems: List<List<Item>>)

}