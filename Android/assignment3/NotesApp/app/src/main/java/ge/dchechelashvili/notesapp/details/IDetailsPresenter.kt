package ge.dchechelashvili.notesapp.details

import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

interface IDetailsPresenter {
    abstract fun saveInfo(noteId: Int, name: String, pinned: Boolean, checkedItems: MutableList<Item>, uncheckedItems: MutableList<Item>)
    abstract fun onInfoSaved()
    abstract fun getNote(id: Int)
    abstract fun onNoteFetched(note: Note)
    abstract fun onNoteItemsFetched(items: List<Item>)
}