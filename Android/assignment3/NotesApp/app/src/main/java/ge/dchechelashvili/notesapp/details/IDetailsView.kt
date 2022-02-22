package ge.dchechelashvili.notesapp.details

import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

interface IDetailsView {
    fun exit()
    fun showNote(note: Note)
    abstract fun showCheckedItems(checkedItems: MutableList<Item>)
    abstract fun showUncheckedItems(uncheckedItems: MutableList<Item>)
}