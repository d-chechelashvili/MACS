package ge.dchechelashvili.notesapp.main

import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

interface IMainPresenter {
    abstract fun onNoteListFetched(notes: List<Note>)
    abstract fun onItemListFetched(items: List<List<Item>>)
}